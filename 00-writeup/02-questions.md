# TODO
1. B. LibreChat 不是  

2. 4.3 3.4  

3. 838 TFLOPS 419 TFLOPS  
    FP8 张量核全速 838 TFLOPS 是 FP16 累加路径；改用 FP32 累加走半速路径（Ada 白皮书先例：4080 的 FP8 389.9/194.9 恰好一半，Blackwell sm_120 实测同样减半）。  
    104.8 是 FP32 CUDA Core 峰值（非张量核），不要混用。  

4. 
    2: **enp189s0f0** : <BROADCAST,MULTICAST,UP,LOWER_UP> mtu 1500 qdisc fq_codel state UP mode DEFAULT group default qlen 1000  
    link/ether 10:c3:ab:50:cb:28 brd ff:ff:ff:ff:ff:ff  
    3: **enp189s0f1**: <NO-CARRIER,BROADCAST,MULTICAST,UP> mtu 1500 qdisc fq_codel state DOWN mode DEFAULT group default qlen 1000
    link/ether 10:c3:ab:50:cb:29 brd ff:ff:ff:ff:ff:ff  
    4: **enp125s0f0**: <NO-CARRIER,BROADCAST,MULTICAST,UP> mtu 1500 qdisc mq state DOWN mode DEFAULT group default qlen 1000
    link/ether 10:c3:ab:50:cb:27 brd ff:ff:ff:ff:ff:ff  
    **firmware-version: 1.8.15.0**  
      
    ***如何找到这些信息:***  
    先申请一个计算节点: `salloc -p kp_run -q kp_run -c 1 -t 00:05:00`  
    然后列出所有的ip link: `ip link show`  
    最后打印详细信息获得硬件版本: `ethtool -i enp189s0f0`  

5.   Entry point address: 0x404440  
    `gcc --version` 拿版本号；`readelf -h $(which gcc) | grep "Entry point address"` 拿入口。  
    gcc 与 g++ 的入口均为 0x404440。

6.  每个NUMA38核心，单NUMA FP64理论算力7539.2 GFLOPS  
    精确算法：每核每周期 128 FLOPs（64 FMA）x 1.55 GHz x 38 核 = 7539.2 GFLOPS。  
    60.3 本身是 304x128x1.55=60.3136 的舍入，精确值应为 60.3136/8=7.5392 TFLOPS。

7.  edaad70ee1670b56aeb090a5ebd038925f60146459b74924e0f8f5cc26eecca4  
    要点：题目要的是 9/1 前该 URL 的历史内容（当时是 PNG，网站 9/1 改版成 SPA 后同一 URL 返回 460 字节 index.html）。  
    9/23 之后直接 wget 得到的 f4e9f91c... 是改版后 HTML 的哈希，判错。  
    正确做法：Wayback Machine 找 9/1 前的快照下载（或 OJ 题目配图），再 sha256sum。Windows 上 `Get-FileHash` 输出大写，提交时转小写。  

8. 每线程64寄存器：50.0%
    每线程65寄存器：37.5%  

9. 0.559 μs 2.71 GiB

    **第一空（decode 时第 4 层 sparse attention 最小时间）:**
    第 4 层 compress_ratio=128（config.json 的 compress_ratios[3]）。
    decode 时 q 为 1 token（64 heads × 512 dim），关注的 KV = window 128 个原生 KV + 524288/128 = 4096 个压缩块 = 4224 个 512 维向量。
    sparse_attn 仅两个 GEMM（kernel.py）:
    - scores = q·kvᵀ: [64,512]×[512,4224] → 2×64×4224×512 = 276,824,064 FLOPs
    - o = softmax(scores)·kv: 同上 → 合计 ≈ 5.537×10⁸ FLOPs
    访存（BF16, 每份数据只进出 HBM 一次）: q 64×512×2 + kv 4224×512×2 + o 64×512×2 + topk_idxs(int32) 4224×4 = 4,473,600 B
    B200 SXM: BF16 峰值 2250 TFLOPS (2.25e15), HBM3e 带宽 8 TB/s (8e12)。
    纯计算 = 5.537e8/2.25e15 = 0.246 μs; 纯访存 = 4.4736e6/8e12 = 0.559 μs。
    计算与访存完全重叠 → t = max = **0.559 μs** (memory-bound, 算术强度仅 ~124 FLOP/B)。
    注意 topk_idxs 是 sparse_attn 的必要输入张量（gather 索引），不可漏计；漏了会得 0.557，判错。

    **第二空（主注意力 KV cache）:**
    每层 kv_cache 形状 [1, window+seq_len/ratio, 512] BF16（model.py）:
    - ratio=0（纯滑窗）: 2 层 × 128×512×2 = 2×131,072 B
    - ratio=4: 21 层 × 131200×512×2 = 21×134,348,800 B
    - ratio=128: 20 层 × 4224×512×2 = 20×4,325,376 B
    合计 2,908,094,464 B = 2.708 GiB ≈ **2.71 GiB**
    （不压缩全量约 21.5 GiB，压缩省了约 8 倍；只计 model.layers[i].attn.kv_cache，indexer/compressor 的缓存不算。）
