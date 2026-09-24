#!/bin/bash
#SBATCH -J miniclash-measure
#SBATCH -p kp_run
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH -t 00:30:00
#SBATCH -o measure-%j.out
#SBATCH -e measure-%j.err

# 在目标机器（Kunpeng，32 核）上测三件事：
#   1. 每个任务的核心秒数 mu —— 决定分数的唯一杠杆
#   2. 线程数 8/16/32 下的并行效率 —— 看掉不掉
#   3. 算例 1 的真实挂钟时间（满分 4.5 s、零分 60 s），跑 3 遍看最坏情况
#
# 用法（必须在 05-miniclash 目录下提交）：
#   cd ~/hellohpc/05-miniclash && sbatch bench/measure.sbatch
#
# 结果全在 stdout，落进 measure-<jobid>.out，可直接整份贴回。

set -u
exec 2>&1                      # 让 bash 内建 time 的输出也进 .out，不跑进 .err
cd "$SLURM_SUBMIT_DIR"

if [ ! -x source_code/run ]; then
	echo "!! 找不到 source_code/run。请在 05-miniclash 目录下提交："
	echo "     cd ~/hellohpc/05-miniclash && sbatch bench/measure.sbatch"
	exit 1
fi

BIN=source_code/run
TIMEFORMAT='   >>> real %3R s | user %3U s | sys %3S s'

echo "############ 0. 机器信息（解释效率曲线要用）############"
echo "hostname   = $(hostname)"
echo "nproc      = $(nproc)"
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-未设置}"
grep -m1 -i 'model name' /proc/cpuinfo \
	|| grep -m1 -i 'cpu part' /proc/cpuinfo \
	|| echo "  (/proc/cpuinfo 里没有型号字段)"
echo "Intel(R) 风格缓存大小: L1d=$(getconf LEVEL1_DCACHE_SIZE) L2=$(getconf LEVEL2_CACHE_SIZE) L3=$(getconf LEVEL3_CACHE_SIZE)"
echo "  ↑ Kunpeng 上 getconf 常返回 0 或不支持，属正常，不用管"
echo "  ↓ 关键：Thread(s) per core 是不是 1（Kunpeng 920 无超线程）"
lscpu 2>/dev/null | grep -i -E 'thread|core|socket|cache|mhz|model name|^cpu\(s\)' || echo "  (没有 lscpu)"

echo
echo "############ 1. 生成算例（目录已存在就跳过，避免覆盖）############"
[ -d case32 ] || python3 utils/gencase.py 32 case32 20260924
[ -d case8  ] || python3 utils/gencase.py  8 case8  20260925
echo "  case32: $(awk 'NF{c++}END{print c+0}' case32/tasks.txt) 个任务"
echo "  case8 : $(awk 'NF{c++}END{print c+0}' case8/tasks.txt) 个任务"

echo
echo "############ 2. 基线：8 任务 / 1 线程（纯串行）############"
echo "  mu = user 时间 / 8。这个数字直接决定能拿几分，最重要的一条。"
export OMP_NUM_THREADS=1
time $BIN case8/tasks.txt

echo
echo "############ 3. 效率曲线：32 任务，线程数 8 / 16 / 32 ############"
echo "  效率 = user / (线程数 x real)，1.0 是满分。看它随线程数掉不掉。"
export OMP_NUM_THREADS=8
time $BIN case32/tasks.txt

export OMP_NUM_THREADS=16
time $BIN case32/tasks.txt

echo "  --- 下面这条就是算例 1 的真实预测：满分 4.5 s，零分 60 s ---"
export OMP_NUM_THREADS=32
time $BIN case32/tasks.txt

echo
echo "############ 4. 最坏情况：32 线程再跑两遍算例 1 ############"
echo "  题面原话：控制最坏情况和控制平均值同样重要。所以要连跑三遍看离散度。"
time $BIN case32/tasks.txt
time $BIN case32/tasks.txt

echo
echo "############ 5. 正确性门禁（换了架构必须重跑，这关不过后面的数字都没意义）############"
python3 utils/verify.py case8/tasks.txt
python3 utils/verify.py case32/tasks.txt

echo
echo "############ 完成。把整个 .out 贴回来即可 ############"
