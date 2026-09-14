#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
static uint64_t next(uint64_t *s) {
  uint64_t x = *s;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  *s = x;
  return x * 2685821657736338717ULL;
}
static int ri(uint64_t *s, int n) { return (int)(next(s) % (uint64_t)n); }
static long trades(int8_t *z, int r, int c, uint64_t *s, int attempts, int *idx,
                   int8_t *old) {
  long moves = 0;
  for (int t = 0; t < attempts; t++) {
    int a = ri(s, r), b = ri(s, r);
    if (a == b)
      continue;
    int d = 0, ones = 0;
    for (int j = 0; j < c; j++) {
      int8_t x = z[a * c + j], y = z[b * c + j];
      if (x >= 0 && y >= 0 && x != y) {
        idx[d] = j;
        old[d] = x;
        ones += x;
        d++;
      }
    }
    if (!ones || ones == d)
      continue;
    for (int j = 0; j < d; j++)
      z[a * c + idx[j]] = 0;
    for (int j = 0; j < ones; j++) {
      int k = j + ri(s, d - j);
      int tmp = idx[j];
      idx[j] = idx[k];
      idx[k] = tmp;
      z[a * c + idx[j]] = 1;
    }
    int changed = 0;
    for (int j = 0; j < d; j++) {
      int col = idx[j];
      if (z[a * c + col] + z[b * c + col] != 1)
        changed = 1;
      z[b * c + col] = 1 - z[a * c + col];
    }
    moves += changed;
  }
  return moves;
}
static void calc(const int8_t *z, int r, int c, int i, int j, double *cov,
                 double *phi) {
  int n = 0, a = 0, b = 0, ab = 0;
  for (int k = 0; k < r; k++) {
    int x = z[k * c + i], y = z[k * c + j];
    if (x < 0 || y < 0)
      continue;
    n++;
    a += x;
    b += y;
    ab += x * y;
  }
  if (n < 12) {
    *cov = NAN;
    *phi = NAN;
    return;
  }
  double p = (double)a / n, q = (double)b / n, v = (double)ab / n - p * q;
  *cov = v * n / (n - 1.0) * 11.0 / 12.0 * 100.0;
  double den = sqrt(p * (1 - p) * q * (1 - q));
  *phi = den > 0 ? v / den : NAN;
}
/* out[p,6]: observed covariance, mean-null covariance, chain difference,
   observed phi, mean-null phi (all 64 draws finite), count of finite null phi
   draws. */
int score(const int8_t *x, int r, int c, const int *ii, const int *jj, int p,
          int encoding, uint64_t seed, double *out, long *info) {
  if (r < 2 || c < 2)
    return -1;
  int8_t *original = malloc((size_t)r * c), *z = malloc((size_t)r * c),
         *old = malloc(c);
  int *idx = malloc(c * sizeof(int));
  int *rows = calloc(r, sizeof(int)), *cols = calloc(c, sizeof(int));
  double *sum = calloc(p, sizeof(double)), *ps = calloc(p, sizeof(double)),
         *chain = calloc(p, sizeof(double));
  int *pn = calloc(p, sizeof(int));
  if (!original || !z || !old || !idx || !rows || !cols || !sum || !ps ||
      !chain || !pn)
    return -2;
  for (int k = 0; k < r * c; k++) {
    int v = x[k];
    original[k] = v < 0 ? -1
                        : (encoding == 0   ? v == 2
                           : encoding == 1 ? v == 1
                                           : v > 0);
    if (original[k] == 1) {
      rows[k / c]++;
      cols[k % c]++;
    }
  }
  for (int k = 0; k < p; k++) {
    calc(original, r, c, ii[k], jj[k], &out[k * 6], &out[k * 6 + 3]);
  }
  info[2] = 1;
  for (int ch = 0; ch < 2; ch++) {
    memcpy(z, original, (size_t)r * c);
    uint64_t state = seed ^ (0x9e3779b97f4a7c15ULL * (ch + 1));
    if (!state)
      state = 1;
    long moves = trades(z, r, c, &state, 50 * r, idx, old);
    for (int t = 0; t < 32; t++) {
      moves += trades(z, r, c, &state, 5 * r, idx, old);
      for (int k = 0; k < p; k++) {
        double cv, ph;
        calc(z, r, c, ii[k], jj[k], &cv, &ph);
        sum[k] += cv;
        chain[k] += (ch == 0 ? 1 : -1) * cv;
        if (isfinite(ph)) {
          ps[k] += ph;
          pn[k]++;
        }
      }
    }
    info[ch] = moves;
    for (int a = 0; a < r; a++) {
      int n = 0;
      for (int b = 0; b < c; b++) {
        int k = a * c + b;
        n += z[k] == 1;
        if ((z[k] < 0) != (original[k] < 0))
          info[2] = 0;
      }
      if (n != rows[a])
        info[2] = 0;
    }
    for (int b = 0; b < c; b++) {
      int n = 0;
      for (int a = 0; a < r; a++)
        n += z[a * c + b] == 1;
      if (n != cols[b])
        info[2] = 0;
    }
  }
  for (int k = 0; k < p; k++) {
    out[k * 6 + 1] = sum[k] / 64;
    out[k * 6 + 2] = chain[k] / 32;
    out[k * 6 + 4] = pn[k] == 64 ? ps[k] / 64 : NAN;
    out[k * 6 + 5] = pn[k];
  }
  free(original);
  free(z);
  free(old);
  free(idx);
  free(rows);
  free(cols);
  free(sum);
  free(ps);
  free(chain);
  free(pn);
  return 0;
}
