#include "kernel.c"

/* write a simple test*/
#include <stdio.h>
#include <stdlib.h>

/* --- 简单的伪随机数（线性同余 LCG），固定种子 → 结果可复现 --- */
static unsigned int rng_state = 12345u;
static float rnd(void) {
  rng_state = rng_state * 1103515245u + 12345u;
  return (float)((rng_state >> 8) & 0xFFFFu) / 32768.0f - 1.0f; /* [-1, 1) */
}

int main(void) {
  const int Q = 37, C = 48, D = 16; /* 取奇数 Q，能揪出边界 bug */
  const size_t n_qd = (size_t)Q * D;
  const size_t n_cd = (size_t)C * D;

  float* points = malloc(n_qd * sizeof(float));
  float* centers = malloc(n_cd * sizeof(float));
  float* trig_vec = malloc(n_cd * sizeof(float));
  float* weights = malloc((size_t)C * sizeof(float));
  float* scales = malloc((size_t)C * sizeof(float));
  float* bias = malloc((size_t)C * sizeof(float));
  float* trig_scale = malloc((size_t)C * sizeof(float));
  float* out = malloc((size_t)Q * sizeof(float));
  double* ref = malloc((size_t)Q * sizeof(double));

  if (!points || !centers || !trig_vec || !weights || !scales || !bias ||
      !trig_scale || !out || !ref) {
    fprintf(stderr, "malloc failed\n");
    return 1;
  }

  /* --- 喂入数据 --- */
  for (size_t i = 0; i < n_qd; i++) points[i] = rnd();
  for (size_t i = 0; i < n_cd; i++) centers[i] = rnd();
  for (size_t i = 0; i < n_cd; i++) trig_vec[i] = rnd();
  for (int c = 0; c < C; c++) {
    weights[c] = rnd();
    scales[c] = 0.5f + 0.5f * fabsf(rnd()); /* 保持为正，贴近题意 */
    bias[c] = rnd();
    trig_scale[c] = 2.0f * rnd();
  }

  /* --- 调用被测内核 --- */
  kernel_field(points, centers, weights, scales, bias, trig_scale, trig_vec,
               out, Q, C, D);

  /* --- 用 double 精度算一份"金标准"参考值 --- */
  for (int q = 0; q < Q; q++) {
    const float* x = points + (size_t)q * D;
    double acc = 0.0;
    for (int c = 0; c < C; c++) {
      const float* mu = centers + (size_t)c * D;
      const float* tv = trig_vec + (size_t)c * D;
      double rr = 0.0, pp = 0.0;
      for (int d = 0; d < D; d++) {
        double delta = (double)x[d] - (double)mu[d];
        rr += delta * delta;
        pp += (double)x[d] * (double)tv[d];
      }
      acc += (double)weights[c] * exp(-(double)scales[c] * rr) +
             (double)bias[c] * sin((double)trig_scale[c] * pp);
    }
    ref[q] = acc;
  }

  /* --- 用官方容差公式打分 ---
     |实际 - 参考| ≤ 1e-5 + 1e-4 · |参考|          */
  double worst = 0.0;
  for (int q = 0; q < Q; q++) {
    double tol = 1e-5 + 1e-4 * fabs(ref[q]);
    double err = fabs((double)out[q] - ref[q]);
    double ratio = err / tol;
    if (ratio > worst) worst = ratio;
  }

  printf("worst error / tolerance = %.4f\n", worst);
  printf("%s\n", worst < 1.0 ? "PASS" : "FAIL");

  free(points);
  free(centers);
  free(trig_vec);
  free(weights);
  free(scales);
  free(bias);
  free(trig_scale);
  free(out);
  free(ref);
  return worst < 1.0 ? 0 : 1;
}
