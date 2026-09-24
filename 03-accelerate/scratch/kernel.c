#include <math.h>
#include <stddef.h>

static inline float fast_exp(float x) {
  double t = (double)x;
  double k = rint(t * 1.4426950408889634);
  double r = t - k * 0.69314718055994529;
  // eᵗ = e^(k·ln2 + r) = 2ᵏ * e^r, r is little. So we do taylor to e^r.
  // Taylor:
  double p = 1.394858082e-03; /* Horner，从最高次 c6 开始 */
  p = p * r + 8.375128891e-03;
  p = p * r + 4.166621827e-02;
  p = p * r + 1.666641548e-01;
  p = p * r + 5.000000108e-01;
  p = p * r + 1.000000038e+00;
  p = p * r + 1.000000000e+00;
  return (float)ldexp(p, (int)k);
}

static inline float fast_sin(float x) {
  if (x > 1e7f || x < -1e7f) return sinf(x); /* 保险分支，实际数据用不到 */
  double t = (double)x;
  double k = rint(t * 0.31830988618379069); /* k = round(t / pi) */
  double r = t - k * 3.1415926535897931;    /* r 落到 [-pi/2, pi/2] */
  double u = r * r;                         /* 奇函数：只含 r 的奇次幂 */
  double p = -2.388949828e-08;
  p = p * u + 2.752528712e-06;
  p = p * u + -1.984086155e-04;
  p = p * u + 8.333330978e-03;
  p = p * u + -1.666666662e-01;
  p = p * u + 1.000000000e+00;
  double s = r * p;
  return (float)(((long)k & 1L) ? -s : s); /* sin(t) = (-1)^k * sin(r) */
}

void kernel_field(const float* restrict points, const float* restrict centers,
                  const float* restrict weights, const float* restrict scales,
                  const float* restrict bias, const float* restrict trig_scale,
                  const float* restrict trig_vec, float* restrict out, int Q,
                  int C, int D) {
#pragma omp parallel for schedule(static)
  for (int q = 0; q < Q; q++) {
    const float* x =
        points + (size_t)q * (size_t)D;  // get pointer of q th start
    double acc = 0.0f;

    for (int c = 0; c < C; c++) {
      const float* mu = centers + (size_t)c * (size_t)D;
      const float* tv = trig_vec + (size_t)c * (size_t)D;

      double rr = 0.0f;
      double pp = 0.0f;
      for (int d = 0; d < D; d++) {
        float delta = x[d] - mu[d];
        rr += delta * delta;
        pp += x[d] * tv[d];
      }
      double arg_exp = -(double)scales[c] * rr;
      double arg_sin = (double)trig_scale[c] * pp;
      acc += (double)weights[c] * (double)expf((float)arg_exp) +
             (double)bias[c] * (double)sinf((float)arg_sin);
    }
    out[q] = (float)acc;
  }
}
