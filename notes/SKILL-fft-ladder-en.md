---
name: cpu-hpc-skill
description: fft optimization ladder - twiddles, cross-batch Stockham, radix-4.
---

# fft optimization ladder

S<2 scores 0, S=4 is full marks. Climb one rung at a time.

**Validate each rung with exactly `bash tools/test_candidate.sh`** - it builds, checks correctness and benchmarks in one call, printing `median_time_ms`.

**Rung 1 - twiddle tables.** Precompute one contiguous `cos[]`/`sin[]` table per stage, and a `rev[]` table for bit reversal. (~1.7x)

**Rung 2 - SoA + vectorize.** Split complex into `re[]`/`im[]` so the butterfly auto-vectorizes. (~2.4x, first scoring rung)

**Rung 3 - cross-batch layout + Stockham.** Vectorize across batches, not across butterflies. A call does 16/8 transforms, so lay out `re[i*8+l]` (i=index, l=batch): one row of 8 = 2 NEON vectors, 4 lanes share one twiddle, no shuffles. 8 lanes is the sweet spot (16 is slower, L1). Drop the separate bit-reversal pass - Stockham folds reordering into the stage addresses. Stage s, `h=2^(s-1)`, `j<N/(2h)`, `k<h`: `u=cur[j*h+k]`, `v=cur[j*h+k+N/2]`; `nxt[j*2*h+k]=u+w*v`, `nxt[j*2*h+k+h]=u-w*v`. Indices are rows. (~3.3x)

**Rung 4 - radix-4.** Traps: output offsets 0,2g,g,3g (1,2 swapped); third branch's twiddle is `-i*w1`, not w3. (~3.5x, full marks)

Rules: single-threaded; no runtime allocation (static arrays, no vector/new); no FFT libraries.

If a rung fails correctness twice, write the last passing source back and stop.
