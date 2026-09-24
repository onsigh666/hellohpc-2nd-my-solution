import ctypes
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent   # 仓库根目录

# ── ① 加载共享库 ───────────────────────────────
lib = ctypes.CDLL(str(ROOT / "scratch" / "kernel.so"))

# ── ② 声明参数类型（漏掉这一步 = 段错误）─────────
F32 = ctypes.POINTER(ctypes.c_float)

lib.kernel_field.argtypes = [
    F32,            # points     (Q, D)
    F32,            # centers    (C, D)
    F32,            # weights    (C,)
    F32,            # scales     (C,)
    F32,            # bias       (C,)
    F32,            # trig_scale (C,)
    F32,            # trig_vec   (C, D)
    F32,            # out        (Q,)
    ctypes.c_int,   # Q
    ctypes.c_int,   # C
    ctypes.c_int,   # D
]
lib.kernel_field.restype = None    # C 函数返回 void

# ── ③ 载入真实数据 ─────────────────────────────
case = np.load(ROOT / "data" / "sample.npz")
points = np.ascontiguousarray(case["points"], dtype=np.float32)
centers = np.ascontiguousarray(case["centers"], dtype=np.float32)
weights = np.ascontiguousarray(case["weights"], dtype=np.float32)
scales = np.ascontiguousarray(case["scales"], dtype=np.float32)
bias = np.ascontiguousarray(case["bias"], dtype=np.float32)
trig_scale = np.ascontiguousarray(case["trig_scale"], dtype=np.float32)
trig_vec = np.ascontiguousarray(case["trig_vec"], dtype=np.float32)

Q, D = points.shape
C, _ = centers.shape

out = np.empty(Q, dtype=np.float32)

# ── ④ 调用 C 内核 ──────────────────────────────
lib.kernel_field(
    points.ctypes.data_as(F32),
    centers.ctypes.data_as(F32),
    weights.ctypes.data_as(F32),
    scales.ctypes.data_as(F32),
    bias.ctypes.data_as(F32),
    trig_scale.ctypes.data_as(F32),
    trig_vec.ctypes.data_as(F32),
    out.ctypes.data_as(F32),
    Q, C, D,
)

# ── ⑤ 对照官方参考答案 ─────────────────────────
ref = np.loadtxt(ROOT / "data" / "sample.ref.txt", dtype=np.float64)
tol = 1e-5 + 1e-4 * np.abs(ref)
ratio = np.abs(out.astype(np.float64) - ref) / tol
print("sample: worst error / tolerance = %.4f" % ratio.max())
print("sample: PASS" if ratio.max() < 1.0 else "sample: FAIL")
