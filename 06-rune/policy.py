from __future__ import annotations

from rune_scheduler.api.v1 import (
    DispatchAction,
    DispatchDecision,
    Observation,
    Placement,
    PublicScenarioModel,
    WorkerState,
)

# CPU 超载余量（只作用于"凑批失败"后的 strict 全局选择）。
# 8 种子实测：0.0（严格不超载）最优；放宽回填（136.3/270）反而更差。
CPU_HEADROOM: float = 0.0

# 凑批轮允许的 CPU 超载上限（cpu_capacity 的倍数）。
# 真均值口径实测（8 种子扫 1.0~∞，64 种子全 18 复核）：
# 1.25 最优 —— 恰好等于"再塞一个任务"（任务 CPU 多为容量的 1/4），
# 蹭到 setup 共享又不真正超载；2.0 允许塞太多导致挤爆快机器。
# 64 种子官方口径：2.0 得 209.3，1.25 得 228.5 / 337。
PACK_LIMIT: float = 1.25


class Policy:
    """贪心 + 三个关键规则（全部经 64 种子官方口径实测验证）：

    1. strict CPU：普通选择绝不让机器 CPU 超载（超载会让整机任务按比例变慢）；
    2. type-packing：同类型任务优先挤到本决策点已放同类型的机器上
       凑批共享 setup；超载上限 PACK_LIMIT=1.25（约等于"再塞一个任务"）；
    3. 稀缺优先：能接机器少的类型先挑机器，稳定排序保持到达顺序。
       64 种子官方口径：L=2.0 无排序 209.3 → L=1.25 228.5 → +排序 242.8 / 337。
    """

    def __init__(self, model: PublicScenarioModel) -> None:
        self._model = model
        # 静态信息按 id 建索引，方便查询。构造时只做一次，不占决策时间。
        self._task_types = {t.id: t for t in model.task_types}
        self._workers = {w.id: w for w in model.workers}
        self._max_batch = model.safety_limits.max_batch_size
        # 每种任务类型有几台机器能接（capabilities 满足数），越少越稀缺。
        # 只算静态能力集合，构造时一次算好，不占决策时间。
        self._cap_count = {}
        for t in model.task_types:
            caps = set(t.required_capabilities)
            self._cap_count[t.id] = sum(
                1 for w in model.workers if caps.issubset(w.capabilities)
            )

    def _duration_ratio(self, worker_id: str, task_type: str) -> float:
        """查这台机器跑这种任务类型的时长倍率，越小越快（1/4 = 4 倍速）。"""
        for entry in self._workers[worker_id].duration_multipliers:
            if entry.task_type == task_type:
                return entry.ratio.numerator / entry.ratio.denominator
        return 1.0  # 没写倍率就是基准速度

    def choose_placements(self, observation: Observation) -> DispatchDecision:
        # --- 第一步：挑出"现在能接活"的机器，记录它还剩多少内存 ---
        free_memory: dict[str, int] = {}
        states: dict[str, WorkerState] = {}
        for state in observation.workers:
            if not state.available:
                continue  # 停机中的机器绝对不能派
            spec = self._workers[state.worker_id]
            free_memory[state.worker_id] = spec.memory_capacity - state.memory_reserved
            states[state.worker_id] = state

        placements: list[Placement] = []
        used_memory: dict[str, int] = {}  # 本决策点内，每台机器已经被我占掉多少内存
        used_cpu: dict[str, int] = {}  # 本决策点内，每台机器已经被我占掉多少 CPU
        placed_types: dict[str, set[str]] = {}  # worker -> 本决策点已放的任务类型

        # --- 第二步：逐个就绪任务，找一台最合适的机器 ---
        # type-packing：先尝试"本决策点已放同类型的机器"（无视 CPU 上限，
        # 凑批共享 setup）；凑批失败再退回 strict 全局逻辑。
        # 稀缺优先：能接的机器越少的类型越先挑机器（瓶颈类型先占位）。
        # 注意：只按 _cap_count 稳定排序，同稀缺度保持到达顺序；
        # 不能用类型名做 tiebreak——会打乱无稀缺差异的 profile（实测 18 号暴跌）。
        ready = sorted(
            observation.ready_tasks,
            key=lambda t: self._cap_count[t.task_type],
        )
        for task in ready:
            if len(placements) >= self._max_batch:
                break  # 框架硬上限：一次最多提交 max_batch_size 个

            spec_type = self._task_types[task.task_type]

            best_worker: str | None = None
            best_ratio = float("inf")

            # 第一轮：凑批（允许 CPU 超载到 PACK_LIMIT 倍，超过就退回 strict 分流）
            for worker_id in free_memory:
                if task.task_type not in placed_types.get(worker_id, ()):
                    continue
                worker = self._workers[worker_id]
                if not set(spec_type.required_capabilities).issubset(
                    worker.capabilities
                ):
                    continue
                remaining = (
                    free_memory[worker_id] - used_memory.get(worker_id, 0)
                )
                if remaining < spec_type.memory_reservation:
                    continue
                cpu_budget = worker.cpu_capacity * PACK_LIMIT
                cpu_left = (
                    cpu_budget - states[worker_id].cpu_demand
                    - used_cpu.get(worker_id, 0)
                )
                if cpu_left < spec_type.cpu_demand:
                    continue
                ratio = self._duration_ratio(worker_id, spec_type.id)
                if ratio < best_ratio:
                    best_ratio = ratio
                    best_worker = worker_id

            # 第二轮：凑批失败，strict 全局找（含 CPU 上限）
            if best_worker is None:
                for worker_id in free_memory:
                    worker = self._workers[worker_id]
                    if not set(spec_type.required_capabilities).issubset(
                        worker.capabilities
                    ):
                        continue
                    remaining = (
                        free_memory[worker_id] - used_memory.get(worker_id, 0)
                    )
                    if remaining < spec_type.memory_reservation:
                        continue
                    cpu_budget = worker.cpu_capacity * (1.0 + CPU_HEADROOM)
                    cpu_left = (
                        cpu_budget - states[worker_id].cpu_demand
                        - used_cpu.get(worker_id, 0)
                    )
                    if cpu_left < spec_type.cpu_demand:
                        continue
                    ratio = self._duration_ratio(worker_id, spec_type.id)
                    if ratio < best_ratio:
                        best_ratio = ratio
                        best_worker = worker_id

            if best_worker is not None:
                placements.append(Placement(task.task_id, best_worker))
                used_memory[best_worker] = (
                    used_memory.get(best_worker, 0) + spec_type.memory_reservation
                )
                used_cpu[best_worker] = (
                    used_cpu.get(best_worker, 0) + spec_type.cpu_demand
                )
                placed_types.setdefault(best_worker, set()).add(task.task_type)

        # --- 第三步：有活就派，一个都放不下就等到下一个事件 ---
        if placements:
            return DispatchDecision(DispatchAction.DISPATCH, tuple(placements))
        return DispatchDecision(DispatchAction.DEFER, ())
