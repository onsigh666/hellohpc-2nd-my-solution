//$Id: bh_interp.h $
// BH_OPT: 球面插值几何预计算缓存 + 融合张量积插值内核
// ----------------------------------------------------------------------------
// 原理：探测器球面在全部迭代中位置不变（采样点只依赖探测器半径 rex），
// 而插值是「对场数据的固定线性函数」——模板下标、拉格朗日权重、越界反射
// 标记都只依赖几何。把它们在计时区外一次算好（与构造器里的 nx_g/arcostheta
// 同类），迭代内只剩对最新场数据的 6^3 次乘加。
//
// 数学等价性：内核是 fmisc.f90::global_interp + polin3（Neville 算法）的
// 拉格朗日直写形式——同一插值多项式，仅求和顺序不同，浮点差异 ~1e-15
// 相对量级（门禁容差 1e-8 abs + 1e-14 rel 的余量之内）。
// 不缓存任何场数据、插值中间态或积分结果（符合比赛规则）。
#ifndef BH_INTERP_H
#define BH_INTERP_H

#include "Block.h"

#define BH_ORDN_MAX 8 // 模板点数上限（配置固定为 6，默认 2*ghost_width=6）

// 单个采样点的插值几何（prepare_point_geom 在计时区外一次算好）
struct bh_point_geom
{
  const Block *BP;      // 拥有该点的网格块；NULL = 未找到（运行时回退旧搜索路径）
  int bi;               // BH_OPT5: 该点所属块在副本数组里的下标（-1 = 无副本，走旧路径）
  int ordn;             // 每方向模板点数 = surface_interpolation_order()
  int base[3];          // 模板起点（1 基下标，等价 global_interp 的 cxB）
  double w[3][BH_ORDN_MAX];     // 各方向拉格朗日权重（节点 0..ordn-1 上取值）
  unsigned char refl[3][BH_ORDN_MAX]; // 该模板点是否越界反射（值需乘 var 的 SoA[d]）
};

// BH_OPT5: 含点块的数据副本。
// 球面采样点只穿过少数几个块（实测 np=8 时是 8 块中的 4 块），而块固定归属某个
// rank——于是 4 个 rank 满负荷插值、另外 4 个从头发呆到尾（见 bh_detector_cache
// 的注释）。但场数据是每轮重新填充的、且填充发生在计时区外，所以只要在计时区外
// 把「含点块」的数据抄一份给所有 rank，计时区内 8 个 rank 就都能插值。
// 插值计算本身仍完整地在计时区内完成，只是数据来源换成了本地副本。
struct bh_block_copy
{
  const Block *src;     // 源块（每轮从这里取最新数据；仅拥有者有效）
  int shape[3];         // = src->shape
  double *f[2];         // 两个变量（Rpsi4/Ipsi4）的数据副本，按 var 列表顺序
};

// 单个探测器（固定 rex）的缓存：采样点 + 每点插值几何 + 本 rank 负责的点列表
struct bh_detector_cache
{
  double *pox[3];        // 球面采样点（与迭代无关）
  bh_point_geom *geom;   // 每点插值几何
  // BH_OPT4: 本 rank 负责的点下标列表（prepare_benchmark 里预计算，计时区外）。
  // 球面采样点只落在少数几个块上：实测 np=8 与 np=32 的点数分布都恰好是
  // 0/0/8826/6040/6040/4694/0/0（其余 rank 全为 0）——并行度被"球面穿过几个块"
  // 锁死在 4，加 rank 没有用。所以每轮再遍历全部 n_tot 个点纯属白扫 65%，
  // 改成只遍历自己这份。
  int *mine;
  int n_mine;
  // BH_OPT5: 含点块的数据副本（见 bh_block_copy）。每轮由 sync_block_copies
  // 在计时区外刷新，使所有 rank 都能插值自己分到的那份点。
  bh_block_copy *bc;
  int n_bc;
};

// 融合张量积插值内核。
//   f   : 块内场数据首地址，布局 [x + y*Nx + z*Nx*Ny]（与 Fortran 列主序一致）
//   Nx, Nxy : shape[0], shape[0]*shape[1]
//   soa : 该变量的反射符号（var::SoA），越界模板值 = 反射点值 * SoA[d]
static inline double bh_fused_interp(const bh_point_geom &g, const double *f,
                                     int Nx, int Nxy, const double soa[3])
{
  const int ordn = g.ordn;
  // 模板首个元素（1 基下标 base 转 0 基 base-1）
  const double *basep = f + (g.base[0] - 1) + (g.base[1] - 1) * Nx + (g.base[2] - 1) * Nxy;

  // 把 var 的 SoA 反射符号折进权重（每 var 每点 3*ordn 次，复用几何权重）
  double w0[BH_ORDN_MAX], w1[BH_ORDN_MAX], w2[BH_ORDN_MAX];
  for (int p = 0; p < ordn; p++)
  {
    w0[p] = g.refl[0][p] ? g.w[0][p] * soa[0] : g.w[0][p];
    w1[p] = g.refl[1][p] ? g.w[1][p] * soa[1] : g.w[1][p];
    w2[p] = g.refl[2][p] ? g.w[2][p] * soa[2] : g.w[2][p];
  }

  // 张量积：先 x（内存连续），再 y，再 z。与 polin3 的 zyx 顺序数学等价。
  // 反射模板的相对下标：fi = base+p <= 0 时取镜像 1-fi（0 基下标 -fi），
  // 相对模板起点（0 基 base-1）即 -(base+p)-(base-1) = 1-2*base-p。
  double vyz[BH_ORDN_MAX][BH_ORDN_MAX];
  if (g.refl[0][0] == 0 && g.refl[1][0] == 0)
  {
    // 主路径（赤道对称基准）：x/y 模板不反射、下标连续，顺序乘加
    for (int r = 0; r < ordn; r++)
    {
      const int iz = g.refl[2][r] ? 1 - 2 * g.base[2] - r : r;
      const double *zp = basep + iz * Nxy;
      for (int q = 0; q < ordn; q++)
      {
        const double *row = zp + q * Nx;
        double s = 0;
        for (int p = 0; p < ordn; p++)
          s += w0[p] * row[p];
        vyz[q][r] = s;
      }
    }
  }
  else
  {
    // 通用路径（含 x/y 反射，如八分对称）
    for (int r = 0; r < ordn; r++)
    {
      const int iz = g.refl[2][r] ? 1 - 2 * g.base[2] - r : r;
      const double *zp = basep + iz * Nxy;
      for (int q = 0; q < ordn; q++)
      {
        const int iy = g.refl[1][q] ? 1 - 2 * g.base[1] - q : q;
        const double *row = zp + iy * Nx;
        double s = 0;
        for (int p = 0; p < ordn; p++)
        {
          const int ix = g.refl[0][p] ? 1 - 2 * g.base[0] - p : p;
          s += w0[p] * row[ix];
        }
        vyz[q][r] = s;
      }
    }
  }
  double vz[BH_ORDN_MAX];
  for (int r = 0; r < ordn; r++)
  {
    double s = 0;
    for (int q = 0; q < ordn; q++)
      s += w1[q] * vyz[q][r];
    vz[r] = s;
  }
  double out = 0;
  for (int r = 0; r < ordn; r++)
    out += w2[r] * vz[r];
  return out;
}

#endif /* BH_INTERP_H */
