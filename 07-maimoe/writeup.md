# 07-maimoe 优化总结

## 题目与评分

- 任务：MaiMoe 谱面编码。host 进程（maimoe）把谱面（chart）经共享内存 + futex 派发给多个 kernel worker 进程（maimoe_kernel），每个谱面做 SHA-256 哈希 + 编码后写回输出。
- 约束：只允许改 solution/Kernel.cpp（kernel 侧入口）和 solution/KstroParam.toml（运行参数），host 与引擎源码只读。
- 评分：score = min(max((base−t)/(base−full), 0), 1)² × weight。四个算例的 base / full / 权重：Sample1 400/100/10、Sample2 2333/160/20、Sample3 2333/165/20、Sample4 11451/650/50。平方公式意味着越接近 full 涨分越快；当前水平下 Sample4 每 +1 分还需再快约 150ms，最贵。

## 优化 1：配置调参（toml，零风险高回报）

- **cpu_workers = 8**（上限）：8 个 kernel worker 并行编码，最大的一笔（本地约 0.2 → 53.4 分）。
- **pipeline_mib = 64**（原 1）：host 有「在途字节预算」（CreditBudget）——所有计算线程手里未完成的谱面字节总和不能超过该预算。默认 1MB 对大谱面算例等于把 8 线程的并行度掐死。提到 64MB 解除限流，本机 Sample4 −6.2%，集群 −3.9%。
- **schedule_policy = size_desc**（原 manifest）：谱面按字节数降序调度，大谱面先开工、小谱面最后塞缝，减少 8 个 worker 收尾阶段的尾部等待（tail effect）。集群三轮对比：

| 配置（集群 hellohpc test） | 总分 | 备注 |
|---|---|---|
| manifest（基线） | 64.74 | |
| work_lpt | 64.86 | 按 work_units 排序收效甚微 |
| size_desc | **64.99** | 采用：Sample4 3737.8 → 3702.6 ms |
| size_desc + write_batch=8 | 64.79 | 负优化，回滚 |

## 优化 2：Kernel.cpp 直通化

原 Kernel.cpp 在调用引擎前做了大量防御性/冗余工作：verify_integrity（对输入谱面重算两遍 SHA-256）、defensive_scan、逐字符文本拷贝、输出 memset 等。引擎本身是确定性的、maimoe_check 才是正确性权威，这些全是白干。全部删掉，直接调用 engine::process_chart_into 直通写入，还省掉逐谱面堆分配。字节等价、全部 checker 通过，本地约 53.4 → 61 分。

## 走过的弯路（教训）

- **环境差异**：WSL 本地创建 worker 进程约 10ms/个，Sample1 本地 195ms；集群约 1ms/个，Sample1 87.5ms 直接满分。涉及进程创建的结论不能以本地为准。
- **work_lpt 预测翻车**：本机模拟预测 Sample3 能省 9%，集群实测只省 0.2%——work_units 排序不是这个数据的好代理。
- **write_batch=8 直觉翻车**：以为攒批写能省系统调用，实测更慢——输出全走页缓存（page cache），syscall 不是瓶颈。
- **SHA-NI 半途而废**：ARMv8 硬件 SHA-256 指令吞吐实测是软件实现的 5.5×（875 vs 160 MB/s），但正确性自检 4 轮未过，暂时放弃。

## 最终成绩

- 集群 hellohpc test：64.99 / 100
- OJ 最终评测：**72 / 100**
