import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import solver

case = np.load(ROOT / "data" / "public-large-b.npz")
args = {k: np.ascontiguousarray(case[k], dtype=np.float32) for k in case}
args["points"] = args["points"][:16384]   # 大题的前 16384 个点

t = time.perf_counter()
out = solver.compute_field(**args)
print(f"elapsed: {time.perf_counter() - t:.3f} s")