from __future__ import annotations

import ctypes
import subprocess
import sys
import tempfile
from pathlib import Path

import math
from typing import Any

import numpy as np

# -------内嵌C内核源码----------
_C_SOURCE =r'''

#include <math.h>
#include <stddef.h>

void kernel_field(const float* restrict points,    
                  const float* restrict centers,    
                  const float* restrict weights,    
                  const float* restrict scales,    
                  const float* restrict bias,       
                  const float* restrict trig_scale,
                  const float* restrict trig_vec,   
                  float* restrict out,              
                  int Q, int C, int D) {
#pragma omp parallel for schedule(static)
  for (int q = 0; q < Q; q++) {
    const float* x =
        points + (size_t)q * (size_t)D;  // get pointer of q th start
    double acc = 0.0f;
#pragma omp simd reduction(+:acc)
    for (int c = 0; c < C; c++) {
      const float* mu = centers + (size_t)c * (size_t)D;
      const float* tv = trig_vec + (size_t)c * (size_t)D;

      double rr = 0.0f; 
      double pp = 0.0f; 
#pragma GCC unroll 8
      for (int d = 0; d < D; d++) {
        float delta = x[d] - mu[d];
        rr += delta * delta;
        pp += x[d] * tv[d];
      }
      double arg_exp=-(double)scales[c]*rr;
      double arg_sin=(double)trig_scale[c]*pp;
      acc +=(double) weights[c] * (double)expf((float)arg_exp) +
           (double)bias[c] *(double)sinf((float)arg_sin);
    }
    out[q] =(float)acc;
  }
}

'''

_F32=ctypes.POINTER(ctypes.c_float)

def _build_kernel():
    with tempfile.TemporaryDirectory(prefix="kernel-field-") as tmp:
        src=Path(tmp)/"kernel.c"
        so = Path(tmp) / "kernel.so"
        src.write_text(_C_SOURCE, encoding="utf-8")
        subprocess.run(
            ["gcc", "-O2", "-shared", "-fPIC","-fopenmp","-funroll-loops",
             "-o", str(so), str(src), "-lm"],
            check=True, capture_output=True,
        )
        lib = ctypes.CDLL(str(so))
    lib.kernel_field.argtypes = [
        _F32, _F32, _F32, _F32, _F32, _F32, _F32, _F32,
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ]
    lib.kernel_field.restype = None
    return lib

try:
    _lib = _build_kernel()
except Exception as exc:          # 任何失败都降级，绝不让 import 崩溃
    print(f"kernel-field: C kernel build failed ({exc}); using numpy fallback",
          file=sys.stderr)
    _lib = None








def compute_field(
    points: Any,
    centers: Any,
    weights: Any,
    scales: Any,
    bias: Any,
    trig_scale: Any,
    trig_vec: Any,
) -> np.ndarray:
    """Contestant entry point with a correct, intentionally slow baseline."""
    q_count, dimension = _shape_2d(points, "points")
    c_count, center_dimension = _shape_2d(centers, "centers")
    trig_rows, trig_dimension = _shape_2d(trig_vec, "trig_vec")
    _shape_1d(weights, "weights", expected=c_count)
    _shape_1d(scales, "scales", expected=c_count)
    _shape_1d(bias, "bias", expected=c_count)
    _shape_1d(trig_scale, "trig_scale", expected=c_count)
    if center_dimension != dimension:
        raise ValueError(f"centers has dim {center_dimension}, expected {dimension}")
    if trig_rows != c_count or trig_dimension != dimension:
        raise ValueError(
            f"trig_vec must have shape ({c_count}, {dimension}), "
            f"got ({trig_rows}, {trig_dimension})"
        )
    if _lib is None:
        return _compute_fallback(points,centers,weights,scales,bias,trig_scale,trig_vec)
    # judge 保证float32 且C连续
    points=np.ascontiguousarray(points,dtype=np.float32)
    centers = np.ascontiguousarray(centers, dtype=np.float32)
    weights = np.ascontiguousarray(weights, dtype=np.float32)
    scales = np.ascontiguousarray(scales, dtype=np.float32)
    bias = np.ascontiguousarray(bias, dtype=np.float32)
    trig_scale = np.ascontiguousarray(trig_scale, dtype=np.float32)
    trig_vec = np.ascontiguousarray(trig_vec, dtype=np.float32)

    output = np.empty(q_count, dtype=np.float32)
    
    _lib.kernel_field( 
        points.ctypes.data_as(_F32),
        centers.ctypes.data_as(_F32),
        weights.ctypes.data_as(_F32),
        scales.ctypes.data_as(_F32),
        bias.ctypes.data_as(_F32),
        trig_scale.ctypes.data_as(_F32),
        trig_vec.ctypes.data_as(_F32),
        output.ctypes.data_as(_F32),
        q_count, c_count, dimension,)
    # Original 3 layers of loops version:
    # for q_index in range(q_count):
    #     point = points[q_index]
    #     accumulator = 0.0
    #     for c_index in range(c_count):
    #         center = centers[c_index]
    #         trig_row = trig_vec[c_index]
    #         square_distance = 0.0
    #         dot_product = 0.0
    #         for d_index in range(dimension):
    #             point_value = float(point[d_index])
    #             delta = point_value - float(center[d_index])
    #             square_distance += delta * delta
    #             dot_product += point_value * float(trig_row[d_index])
    #         accumulator += float(weights[c_index]) * math.exp(
    #             -float(scales[c_index]) * square_distance
    #         )
    #         accumulator += float(bias[c_index]) * math.sin(
    #             float(trig_scale[c_index]) * dot_product
    #         )
    #     output[q_index] = accumulator
    return output

def _compute_fallback(points, centers, weights, scales, bias, trig_scale,
                      trig_vec):
    """纯 numpy 分块实现：只在 C 编译失败时使用。"""
    points = np.ascontiguousarray(points, dtype=np.float32)
    centers = np.ascontiguousarray(centers, dtype=np.float32)
    weights = np.ascontiguousarray(weights, dtype=np.float32)
    scales = np.ascontiguousarray(scales, dtype=np.float32)
    bias = np.ascontiguousarray(bias, dtype=np.float32)
    trig_scale = np.ascontiguousarray(trig_scale, dtype=np.float32)
    trig_vec = np.ascontiguousarray(trig_vec, dtype=np.float32)

    q_count = points.shape[0]
    output = np.empty(q_count, dtype=np.float32)
    for start in range(0, q_count, 256):
        block = points[start:start + 256]                 # (b, D)
        delta = block[:, None, :] - centers[None, :, :]   # (b, C, D)
        rr = (delta * delta).sum(axis=2)                  # 平方距离
        pp = np.einsum("bd,cd->bc", block, trig_vec)      # 点积
        acc = weights * np.exp(-scales * rr) + bias * np.sin(trig_scale * pp)
        output[start:start + 256] = acc.sum(axis=1)       # 跨中心归约
    return output

def _shape_2d(value: Any, name: str) -> tuple[int, int]:
    shape = getattr(value, "shape", None)
    if shape is not None:
        if len(shape) != 2:
            raise ValueError(f"{name} must be 2D, got shape {shape}")
        return int(shape[0]), int(shape[1])
    rows = len(value)
    columns = len(value[0]) if rows else 0
    if any(len(row) != columns for row in value):
        raise ValueError(f"{name} must be rectangular")
    return rows, columns


def _shape_1d(value: Any, name: str, expected: int | None = None) -> int:
    shape = getattr(value, "shape", None)
    if shape is not None:
        if len(shape) != 1:
            raise ValueError(f"{name} must be 1D, got shape {shape}")
        size = int(shape[0])
    else:
        size = len(value)
    if expected is not None and size != expected:
        raise ValueError(f"{name} must have length {expected}, got {size}")
    return size
