import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import solver

print("kernel loaded:", solver._lib is not None)


def check(name, out, ref):
    ref = np.asarray(ref, dtype=np.float64)
    tol = 1e-5 + 1e-4 * np.abs(ref)
    ratio = np.abs(out.astype(np.float64) - ref) / tol
    print(f"{name}: worst error / tolerance = {ratio.max():.4f}  "
          + ("PASS" if ratio.max() < 1.0 else "FAIL"))


# ① sample（完整，走 C 路径）
case = np.load(ROOT / "data" / "sample.npz")
args = {k: np.ascontiguousarray(case[k], dtype=np.float32) for k in case}
out = solver.compute_field(**args)
ref = np.loadtxt(ROOT / "data" / "sample.ref.txt")
check("sample (C)", out, ref)

# ② large-b 抽查（前 2048 点 × 全部 10240 中心，覆盖 D=32 路径）
case_b = np.load(ROOT / "data" / "public-large-b.npz")
args_b = {k: np.ascontiguousarray(case_b[k], dtype=np.float32) for k in case_b}
args_b["points"] = args_b["points"][:2048]
out_b = solver.compute_field(**args_b)
ref_b = np.loadtxt(ROOT / "data" / "public-large-b.ref.txt")[:2048]
check("large-b subset (C)", out_b, ref_b)

# ③ 降级路径也必须测：强制关掉 C，走 numpy fallback
solver._lib = None
out_f = solver.compute_field(**args)
check("sample (fallback)", out_f, ref)
