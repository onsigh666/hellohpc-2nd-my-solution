#include <iostream>

#include <cstdint>

#include "main.hpp"

thread_local uint32 seed32_1, seed32_2;

// 叫停开关：由调用方（批量模式的调度器）设置成本线程正在做的那个任务的
// "已解决"标志。为空表示不叫停（单碰撞 CLI 路径就是这种情况）。
thread_local const std::atomic<bool>* g_stop = nullptr;

namespace {

// 三张差分掩码表只取决于下标 k，与 IV、block 都无关，是纯常量。
// 原来在 find_block0 里每次调用现分配现填(共 256 KB + 64K 次填表)，
// 现在只在程序启动时构造一次。
struct mask_tables {
  uint32 q4[1 << 4];
  uint32 q9q10[1 << 3];
  // SIMD 版 counter4 循环按 4 个一组做 16 字节加载，表头必须 16 字节对齐。
  alignas(16) uint32 q9[1 << 16];

  mask_tables() {
    for (unsigned k = 0; k < (1 << 4); ++k)
      q4[k] = ((k << 2) ^ (k << 26)) & 0x38000004;
    for (unsigned k = 0; k < (1 << 3); ++k)
      q9q10[k] = ((k << 13) ^ (k << 4)) & 0x2060;
    for (unsigned k = 0; k < (1 << 16); ++k)
      q9[k] = ((k << 1) ^ (k << 2) ^ (k << 5) ^ (k << 7) ^ (k << 8) ^
               (k << 10) ^ (k << 11) ^ (k << 13)) &
              0x0eb94f16;
  }
};

// ---------------- counter4 循环的 SIMD 工具（GCC/Clang 向量扩展） ----------------
// u32x4 = 4 个 uint32 的向量：x86_64 映射到 SSE2、AArch64 映射到 NEON，
// 同一份源码两端都能编，本地对拍即验真身（服务器用 NEON，寄存器 32 个）。
typedef uint32 u32x4 __attribute__((vector_size(16)));
typedef int32_t s32x4 __attribute__((vector_size(16)));

inline u32x4 vdup(uint32 x) { return u32x4{x, x, x, x}; }
// 循环左/右移。NEON 没有 32 位 ror 指令，向量版 = 2 移 1 或，这是向量每步
// 比标量(1 条 ROR)贵的主要来源，其余运算向量/标量指令数相当。
inline u32x4 vrotl(u32x4 x, unsigned n) { return (x << n) | (x >> (32 - n)); }
inline u32x4 vrotr(u32x4 x, unsigned n) { return (x >> n) | (x << (32 - n)); }
// 每 lane 算术右移 31 位：该 lane 的 bit31 为 1 则得全 1(死)，否则全 0(活)。
inline u32x4 vshra31(u32x4 x) { return (u32x4)((s32x4)x >> 31); }
// alive 每 lane 全 1=活、全 0=死。用 64 位视角横向或，只需 4 条标量指令。
inline bool all_dead(u32x4 alive) {
  union {
    u32x4 v;
    uint64_t u[2];
  } tmp = {alive};
  return !(tmp.u[0] | tmp.u[1]);
}
inline u32x4 vFF(u32x4 b, u32x4 c, u32x4 d) { return d ^ (b & (c ^ d)); }
inline u32x4 vGG(u32x4 b, u32x4 c, u32x4 d) { return c ^ (d & (b ^ c)); }
inline u32x4 vHH(u32x4 b, u32x4 c, u32x4 d) { return b ^ c ^ d; }
inline u32x4 vII(u32x4 b, u32x4 c, u32x4 d) { return c ^ (b | ~d); }

#define VSTEP(f, a, b, c, d, m, ac, rc) \
  (a = vrotl((a) + f(b, c, d) + (m) + vdup(ac), rc) + (b))

// 文件作用域静态对象：main 之前构造一次，之后全程只读。
// 只读 => 所有搜索线程共享同一份，无数据竞争，也不需要加锁。
const mask_tables g_masks;

}  // namespace

bool find_block0(uint32 block[], const uint32 IV[]) {
  uint32 Q[68] = {IV[0], IV[3], IV[2], IV[1]};

  while (true) {
    // 每轮外层循环查一次叫停开关。距上次检查最多隔一个 2^16 内层枚举
    // (几毫秒)，所以叫停延迟可以忽略。relaxed 足够：只需"最终能看到"。
    if (g_stop && g_stop->load(std::memory_order_relaxed)) return false;

    Q[Qoff + 1] = xrng64();
    Q[Qoff + 3] = (xrng64() & 0xfe87bc3f) | 0x017841c0;
    Q[Qoff + 4] =
        (xrng64() & 0x44000033) | 0x000002c0 | (Q[Qoff + 3] & 0x0287bc00);
    Q[Qoff + 5] = 0x41ffffc8 | (Q[Qoff + 4] & 0x04000033);
    Q[Qoff + 6] = 0xb84b82d6;
    Q[Qoff + 7] = (xrng64() & 0x68000084) | 0x02401b43;
    Q[Qoff + 8] =
        (xrng64() & 0x2b8f6e04) | 0x005090d3 | (~Q[Qoff + 7] & 0x40000000);
    Q[Qoff + 9] =
        0x20040068 | (Q[Qoff + 8] & 0x00020000) | (~Q[Qoff + 8] & 0x40000000);
    Q[Qoff + 10] = (xrng64() & 0x40000000) | 0x1040b089;
    Q[Qoff + 11] =
        (xrng64() & 0x10408008) | 0x0fbb7f16 | (~Q[Qoff + 10] & 0x40000000);
    Q[Qoff + 12] =
        (xrng64() & 0x1ed9df7f) | 0x00022080 | (~Q[Qoff + 11] & 0x40200000);
    Q[Qoff + 13] = (xrng64() & 0x5efb4f77) | 0x20049008;
    Q[Qoff + 14] =
        (xrng64() & 0x1fff5f77) | 0x0000a088 | (~Q[Qoff + 13] & 0x40000000);
    Q[Qoff + 15] =
        (xrng64() & 0x5efe7ff7) | 0x80008000 | (~Q[Qoff + 14] & 0x00010000);
    Q[Qoff + 16] =
        (xrng64() & 0x1ffdffff) | 0xa0000000 | (~Q[Qoff + 15] & 0x40020000);

    MD5_REVERSE_STEP(0, 0xd76aa478, 7);
    MD5_REVERSE_STEP(6, 0xa8304613, 17);
    MD5_REVERSE_STEP(7, 0xfd469501, 22);
    MD5_REVERSE_STEP(11, 0x895cd7be, 22);
    MD5_REVERSE_STEP(14, 0xa679438e, 17);
    MD5_REVERSE_STEP(15, 0x49b40821, 22);

    const uint32 tt1 =
        FF(Q[Qoff + 1], Q[Qoff + 0], Q[Qoff - 1]) + Q[Qoff - 2] + 0xe8c7b756;
    const uint32 tt17 = GG(Q[Qoff + 16], Q[Qoff + 15], Q[Qoff + 14]) +
                        Q[Qoff + 13] + 0xf61e2562;
    const uint32 tt18 = Q[Qoff + 14] + 0xc040b340 + block[6];
    const uint32 tt19 = Q[Qoff + 15] + 0x265e5a51 + block[11];
    const uint32 tt20 = Q[Qoff + 16] + 0xe9b6c7aa + block[0];
    const uint32 tt5 = RR(Q[Qoff + 6] - Q[Qoff + 5], 12) -
                       FF(Q[Qoff + 5], Q[Qoff + 4], Q[Qoff + 3]) - 0x4787c62a;

    // change q17 until conditions are met on q18, q19 and q20
    unsigned counter = 0;
    while (counter < (1 << 7)) {
      const uint32 q16 = Q[Qoff + 16];
      uint32 q17 = ((xrng64() & 0x3ffd7ff7) | (q16 & 0xc0008008)) ^ 0x40000000;
      ++counter;

      uint32 q18 = GG(q17, q16, Q[Qoff + 15]) + tt18;
      q18 = RL(q18, 9);
      q18 += q17;
      if (0x00020000 != ((q18 ^ q17) & 0xa0020000)) continue;

      uint32 q19 = GG(q18, q17, q16) + tt19;
      q19 = RL(q19, 14);
      q19 += q18;
      if (0x80000000 != (q19 & 0x80020000)) continue;

      uint32 q20 = GG(q19, q18, q17) + tt20;
      q20 = RL(q20, 20);
      q20 += q19;
      if (0x00040000 != ((q20 ^ q19) & 0x80040000)) continue;

      block[1] = q17 - q16;
      block[1] = RR(block[1], 5);
      block[1] -= tt17;
      uint32 q2 = block[1] + tt1;
      q2 = RL(q2, 12);
      q2 += Q[Qoff + 1];
      block[5] = tt5 - q2;

      Q[Qoff + 2] = q2;
      Q[Qoff + 17] = q17;
      Q[Qoff + 18] = q18;
      Q[Qoff + 19] = q19;
      Q[Qoff + 20] = q20;
      MD5_REVERSE_STEP(2, 0x242070db, 17);

      counter = 0;
      break;
    }
    if (counter != 0) continue;

    const uint32 q4 = Q[Qoff + 4];
    const uint32 q9backup = Q[Qoff + 9];
    const uint32 tt21 = GG(Q[Qoff + 20], Q[Qoff + 19], Q[Qoff + 18]) +
                        Q[Qoff + 17] + 0xd62f105d;

    // iterate over possible changes of q4
    // while keeping all conditions on q1-q20 intact
    // this changes m3, m4, m5 and m7
    unsigned counter2 = 0;
    while (counter2 < (1 << 4)) {
      Q[Qoff + 4] = q4 ^ g_masks.q4[counter2];
      ++counter2;
      MD5_REVERSE_STEP(5, 0x4787c62a, 12);
      uint32 q21 = tt21 + block[5];
      q21 = RL(q21, 5);
      q21 += Q[Qoff + 20];
      if (0 != ((q21 ^ Q[Qoff + 20]) & 0x80020000)) continue;

      Q[Qoff + 21] = q21;
      MD5_REVERSE_STEP(3, 0xc1bdceee, 22);
      MD5_REVERSE_STEP(4, 0xf57c0faf, 7);
      MD5_REVERSE_STEP(7, 0xfd469501, 22);

      const uint32 tt22 = GG(Q[Qoff + 21], Q[Qoff + 20], Q[Qoff + 19]) +
                          Q[Qoff + 18] + 0x02441453;
      const uint32 tt23 = Q[Qoff + 19] + 0xd8a1e681 + block[15];
      const uint32 tt24 = Q[Qoff + 20] + 0xe7d3fbc8 + block[4];

      const uint32 tt9 = Q[Qoff + 6] + 0x8b44f7af;
      const uint32 tt10 = Q[Qoff + 7] + 0xffff5bb1;
      const uint32 tt8 =
          FF(Q[Qoff + 8], Q[Qoff + 7], Q[Qoff + 6]) + Q[Qoff + 5] + 0x698098d8;
      const uint32 tt12 = RR(Q[Qoff + 13] - Q[Qoff + 12], 7) - 0x6b901122;
      const uint32 tt13 = RR(Q[Qoff + 14] - Q[Qoff + 13], 12) -
                          FF(Q[Qoff + 13], Q[Qoff + 12], Q[Qoff + 11]) -
                          0xfd987193;

      // iterate over possible changes of q9 and q10
      // while keeping conditions on q1-q21 intact
      // this changes m8, m9, m10, m12 and m13 (and not m11!)
      // the possible changes of q9 that also do not change m10 are used below
      for (unsigned counter3 = 0; counter3 < (1 << 3);) {
        uint32 q10 = Q[Qoff + 10] ^ (g_masks.q9q10[counter3] & 0x60);
        Q[Qoff + 9] = q9backup ^ (g_masks.q9q10[counter3] & 0x2000);
        ++counter3;
        uint32 m10 = RR(Q[Qoff + 11] - q10, 17);
        m10 -= FF(q10, Q[Qoff + 9], Q[Qoff + 8]) + tt10;

        uint32 aa = Q[Qoff + 21];
        uint32 dd = tt22 + m10;
        dd = RL(dd, 9) + aa;
        if (0x80000000 != (dd & 0x80000000)) continue;

        uint32 bb = Q[Qoff + 20];
        uint32 cc = tt23 + GG(dd, aa, bb);
        if (0 != (cc & 0x20000)) continue;
        cc = RL(cc, 14) + dd;
        if (0 != (cc & 0x80000000)) continue;

        bb = tt24 + GG(cc, dd, aa);
        bb = RL(bb, 20) + cc;
        if (0 == (bb & 0x80000000)) continue;

        block[10] = m10;
        block[13] = tt13 - q10;

        // iterate over possible changes of q9
        // while keeping intact conditions on q1-q24
        // this changes m8, m9 and m12 (but not m10!)
        //
        // SIMD 版：一次枚举 4 个 q9 候选。18 个早退检查改成"存活 lane 掩码"：
        // 死 lane 继续算（垃圾值不看），整组全死才换下一组。理论加速上限 =
        // 步数比：标量每候选平均 ~22 步、向量每组 4 个平均 ~30 步，≈2.9×；
        // 实际被向量旋转(2移1或)与全死检查拖累，以实测为准。
        const u32x4 vQ8 = vdup(Q[Qoff + 8]), vQ7 = vdup(Q[Qoff + 7]),
                    vQ9 = vdup(Q[Qoff + 9]), vq10 = vdup(q10);
        const u32x4 v_aa = vdup(aa), v_bb = vdup(bb), v_cc = vdup(cc),
                    v_dd = vdup(dd);
        const u32x4 vb0 = vdup(block[0]), vb1 = vdup(block[1]),
                    vb2 = vdup(block[2]), vb3 = vdup(block[3]),
                    vb4 = vdup(block[4]), vb5 = vdup(block[5]),
                    vb6 = vdup(block[6]), vb7 = vdup(block[7]),
                    vb10 = vdup(block[10]), vb11 = vdup(block[11]),
                    vb13 = vdup(block[13]), vb14 = vdup(block[14]),
                    vb15 = vdup(block[15]);
        // block[12] = (tt12 - FF(Q12, Q11, q10)) - q9：括号部分不随 q9 变。
        const u32x4 vb12_base = vdup(tt12 - FF(Q[Qoff + 12], Q[Qoff + 11], q10));
        const u32x4 vtt8 = vdup(tt8), vtt9 = vdup(tt9);

        for (unsigned counter4 = 0; counter4 < (1 << 16); counter4 += 4) {
          // 4 个候选的 q9 与随 q9 变的三个块词。g_masks.q9 已按 16 字节对齐，
          // counter4 步长是 4，所以这是对齐的 16 字节加载。
          u32x4 q9 = vQ9 ^ *(const u32x4*)&g_masks.q9[counter4];
          u32x4 b12 = vb12_base - q9;
          u32x4 b8 = vrotr(q9 - vQ8, 7) - vtt8;
          u32x4 b9 = vrotr(vq10 - q9, 12) - vFF(q9, vQ8, vQ7) - vtt9;

          u32x4 a = v_aa, b = v_bb, c = v_cc, d = v_dd;
          u32x4 alive = vdup(0xFFFFFFFFu);
          VSTEP(vGG, a, b, c, d, b9, 0x21e1cde6, 5);
          VSTEP(vGG, d, a, b, c, vb14, 0xc33707d6, 9);
          VSTEP(vGG, c, d, a, b, vb3, 0xf4d50d87, 14);
          VSTEP(vGG, b, c, d, a, b8, 0x455a14ed, 20);
          VSTEP(vGG, a, b, c, d, vb13, 0xa9e3e905, 5);
          VSTEP(vGG, d, a, b, c, vb2, 0xfcefa3f8, 9);
          VSTEP(vGG, c, d, a, b, vb7, 0x676f02d9, 14);
          VSTEP(vGG, b, c, d, a, b12, 0x8d2a4c8a, 20);
          VSTEP(vHH, a, b, c, d, vb5, 0xfffa3942, 4);
          VSTEP(vHH, d, a, b, c, b8, 0x8771f681, 11);

          c += vHH(d, a, b) + vb11 + vdup(0x6d9d6122u);
          // 检查 1：c 的 bit15 必须为 0。左移 16 把 bit15 顶到 bit31，
          // vshra31 得全 1(死)/全 0(活)，&=~ 即杀置位 lane。
          alive &= ~vshra31(c << 16);
          if (all_dead(alive)) continue;
          c = vrotl(c, 16) + d;

          VSTEP(vHH, b, c, d, a, vb14, 0xfde5380c, 23);
          VSTEP(vHH, a, b, c, d, vb1, 0xa4beea44, 4);
          VSTEP(vHH, d, a, b, c, vb4, 0x4bdecfa9, 11);
          VSTEP(vHH, c, d, a, b, vb7, 0xf6bb4b60, 16);
          VSTEP(vHH, b, c, d, a, vb10, 0xbebfbc70, 23);
          VSTEP(vHH, a, b, c, d, vb13, 0x289b7ec6, 4);
          VSTEP(vHH, d, a, b, c, vb0, 0xeaa127fa, 11);
          VSTEP(vHH, c, d, a, b, vb3, 0xd4ef3085, 16);
          VSTEP(vHH, b, c, d, a, vb6, 0x04881d05, 23);
          VSTEP(vHH, a, b, c, d, b9, 0xd9d4d039, 4);
          VSTEP(vHH, d, a, b, c, b12, 0xe6db99e5, 11);
          VSTEP(vHH, c, d, a, b, vb15, 0x1fa27cf8, 16);
          VSTEP(vHH, b, c, d, a, vb2, 0xc4ac5665, 23);
          // 检查 2：(b^d) 的 bit31 必须为 0。
          alive &= ~vshra31(b ^ d);
          if (all_dead(alive)) continue;

          // II 段 15 个检查，死因有两种："死 if 置位"用 ~，"死 if 清零"不用。
          VSTEP(vII, a, b, c, d, vb0, 0xf4292244, 6);
          alive &= ~vshra31(a ^ c); if (all_dead(alive)) continue;
          VSTEP(vII, d, a, b, c, vb7, 0x432aff97, 10);
          alive &= vshra31(b ^ d); if (all_dead(alive)) continue;
          VSTEP(vII, c, d, a, b, vb14, 0xab9423a7, 15);
          alive &= ~vshra31(a ^ c); if (all_dead(alive)) continue;
          VSTEP(vII, b, c, d, a, vb5, 0xfc93a039, 21);
          alive &= ~vshra31(b ^ d); if (all_dead(alive)) continue;
          VSTEP(vII, a, b, c, d, b12, 0x655b59c3, 6);
          alive &= ~vshra31(a ^ c); if (all_dead(alive)) continue;
          VSTEP(vII, d, a, b, c, vb3, 0x8f0ccc92, 10);
          alive &= ~vshra31(b ^ d); if (all_dead(alive)) continue;
          VSTEP(vII, c, d, a, b, vb10, 0xffeff47d, 15);
          alive &= ~vshra31(a ^ c); if (all_dead(alive)) continue;
          VSTEP(vII, b, c, d, a, vb1, 0x85845dd1, 21);
          alive &= ~vshra31(b ^ d); if (all_dead(alive)) continue;
          VSTEP(vII, a, b, c, d, b8, 0x6fa87e4f, 6);
          alive &= ~vshra31(a ^ c); if (all_dead(alive)) continue;
          VSTEP(vII, d, a, b, c, vb15, 0xfe2ce6e0, 10);
          alive &= ~vshra31(b ^ d); if (all_dead(alive)) continue;
          VSTEP(vII, c, d, a, b, vb6, 0xa3014314, 15);
          alive &= ~vshra31(a ^ c); if (all_dead(alive)) continue;
          VSTEP(vII, b, c, d, a, vb13, 0x4e0811a1, 21);
          alive &= vshra31(b ^ d); if (all_dead(alive)) continue;
          VSTEP(vII, a, b, c, d, vb4, 0xf7537e82, 6);
          alive &= ~vshra31(a ^ c); if (all_dead(alive)) continue;
          VSTEP(vII, d, a, b, c, vb11, 0xbd3af235, 10);
          alive &= ~vshra31(b ^ d); if (all_dead(alive)) continue;
          VSTEP(vII, c, d, a, b, vb2, 0x2ad7d2bb, 15);
          alive &= ~vshra31(a ^ c); if (all_dead(alive)) continue;
          VSTEP(vII, b, c, d, a, b9, 0xeb86d391, 21);

          // 40 步 + 18 检查全过的 lane（极稀有）。按 lane 升序走原标量终检，
          // 最先成功的立即返回 —— 与标量版"先找到先返回"完全一致，输出逐字节相同。
          union {
            u32x4 v;
            uint32_t u[4];
          } lb = {b}, lc = {c}, ld = {d}, l8 = {b8}, l9 = {b9}, l12 = {b12},
            lm = {alive};
          for (unsigned lane = 0; lane < 4; ++lane) {
            if (0 == lm.u[lane]) continue;
            block[12] = l12.u[lane];
            block[8] = l8.u[lane];
            block[9] = l9.u[lane];

            uint32 IHV1 = lb.u[lane] + IV[1];
            uint32 IHV2 = lc.u[lane] + IV[2];
            uint32 IHV3 = ld.u[lane] + IV[3];

            bool wang = true;
            if (0x02000000 != ((IHV2 ^ IHV1) & 0x86000000)) wang = false;
            if (0 != ((IHV1 ^ IHV3) & 0x82000000)) wang = false;
            if (0 != (IHV1 & 0x06000020)) wang = false;

            bool stevens = true;
            if (((IHV1 ^ IHV2) >> 31) != 0 || ((IHV1 ^ IHV3) >> 31) != 0)
              stevens = false;
            if ((IHV3 & (1 << 25)) != 0 || (IHV2 & (1 << 25)) != 0 ||
                (IHV1 & (1 << 25)) != 0 || ((IHV2 ^ IHV1) & 1) != 0)
              stevens = false;

            if (!(wang || stevens)) continue;

            std::cout << "." << std::flush;

            uint32 IV1[4], IV2[4];
            for (int t = 0; t < 4; ++t) IV2[t] = IV1[t] = IV[t];

            uint32 block2[16];
            for (int t = 0; t < 16; ++t) block2[t] = block[t];
            block2[4] += 1 << 31;
            block2[11] += 1 << 15;
            block2[14] += 1 << 31;

            md5_compress(IV1, block);
            md5_compress(IV2, block2);
            if ((IV2[0] == IV1[0] + (1 << 31)) &&
                (IV2[1] == IV1[1] + (1 << 31) + (1 << 25)) &&
                (IV2[2] == IV1[2] + (1 << 31) + (1 << 25)) &&
                (IV2[3] == IV1[3] + (1 << 31) + (1 << 25)))
              return true;

            if (IV2[0] != IV1[0] + (1 << 31)) std::cout << "!" << std::flush;
          }
        }
      }
    }
  }
}
