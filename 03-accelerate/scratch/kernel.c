#include <math.h>
#include <stddef.h>

void kernel_field(const float* points, const float* centers,
                  const float* weights, const float* scales, const float* bias,
                  const float* trig_scale, const float* trig_vec, float* out,
                  int Q, int C, int D) {
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
