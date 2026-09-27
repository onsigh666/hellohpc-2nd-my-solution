---
name: cpu-hpc-skill
description: fft 内核优化手册：查表、SoA+SIMD、radix-4 三步改法；不到 2× 得 0 分。动手前先加载。
---

# fft 优化手册

计分只看加速比：快 2× 才开始得分，4× 满分，不到 2× 一律 0 分。前两步只到 1.7×，第 3 步是必做项：

1. **位反转查表**：starter 每元素现算 log_n 次移位；预生成 `rev[]` 表，置换只查表。
2. **旋转因子查表**：starter 蝶形里递推累乘；改成每层连续表 `w[j]=exp(-2πi·j/span)`。
3. **SoA+向量化**：复数拆成 `re[]/im[]` 两条数组，蝶形 j 循环加 `#pragma omp simd`。三步做完约 2.4×，没到 2× 别停。
4. 仍不够：radix-4（两层合并、访存减半）、按 batch 展开喂满向量单元。

硬约束：单线程；不得运行时分配内存（用 static 数组，禁 vector/new）；禁外部 FFT 库。

流程：每改一版 → `bash tools/check.sh`（须 `correctness: passed`）→ `bash tools/bench.sh` 读 `median_time_ms` → `bash tools/vec_report.sh` 确认向量化。
