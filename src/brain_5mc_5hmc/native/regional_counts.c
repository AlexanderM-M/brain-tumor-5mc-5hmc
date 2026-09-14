#include "htslib/khash.h"
#include "htslib/sam.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
KHASH_SET_INIT_STR(readnames)
typedef struct {
  const uint32_t *sites, *offsets, *lengths;
  uint32_t *counts;
  uint64_t *regions, *qc;
  khash_t(readnames) * seen;
  uint32_t nsites, nregions;
} context;
/* Sites are zero-based reference CpG C coordinates, autosomes in numeric order.
 * Site fields: C80,m80,h80,failed80,explicit_paired,reference_opportunity,
 * intact_read_CpG_opportunity,missing_at_intact_read_CpG.
 * Regions: 100 kb, 12 state counts (C,m,h at .70,.80,.90,.95).
 * Canonical residual clamped within one quantization bin, as validated pilot.
 */
void *regional_create(const uint32_t *sites, const uint32_t *offsets,
                      const uint32_t *lengths, uint32_t *counts,
                      uint64_t *regions, uint64_t *qc, uint32_t nsites,
                      uint32_t nregions) {
  context *x = calloc(1, sizeof(context));
  if (!x)
    return NULL;
  x->sites = sites;
  x->offsets = offsets;
  x->lengths = lengths;
  x->counts = counts;
  x->regions = regions;
  x->qc = qc;
  x->nsites = nsites;
  x->nregions = nregions;
  x->seen = kh_init(readnames);
  if (!x->seen) {
    free(x);
    return NULL;
  }
  return x;
}
void regional_destroy(void *ptr) {
  context *x = ptr;
  if (!x)
    return;
  khiter_t k;
  for (k = kh_begin(x->seen); k != kh_end(x->seen); ++k)
    if (kh_exist(x->seen, k))
      free((char *)kh_key(x->seen, k));
  kh_destroy(readnames, x->seen);
  free(x);
}
static uint32_t lower(const uint32_t *a, uint32_t lo, uint32_t hi, int64_t v) {
  while (lo < hi) {
    uint32_t m = lo + (hi - lo) / 2;
    if ((int64_t)a[m] < v)
      lo = m + 1;
    else
      hi = m;
  }
  return lo;
}
int regional_process(void *ptr, const char *path, char *err, size_t errlen) {
  context *x = ptr;
  samFile *f = sam_open(path, "r");
  sam_hdr_t *h = NULL;
  bam1_t *b = NULL;
  hts_base_mod_state *state = NULL;
  int status = -1, ret;
#define FAIL(msg)                                                              \
  do {                                                                         \
    snprintf(err, errlen, "%s", msg);                                          \
    goto finish;                                                               \
  } while (0)
  if (!f)
    FAIL("Cannot open BAM");
  h = sam_hdr_read(f);
  b = bam_init1();
  state = hts_base_mod_state_alloc();
  if (!h || !b || !state)
    FAIL("Allocation/header error");
  uint32_t region_offsets[22];
  region_offsets[0] = 0;
  for (int c = 1; c < 22; c++)
    region_offsets[c] =
        region_offsets[c - 1] + (x->lengths[c - 1] + 99999) / 100000;
  while ((ret = sam_read1(f, h, b)) >= 0) {
    x->qc[0]++;
    uint16_t flag = b->core.flag;
    if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
      x->qc[1]++;
      continue;
    }
    if (flag & BAM_FUNMAP) {
      x->qc[2]++;
      continue;
    }
    if (flag & (BAM_FDUP | BAM_FQCFAIL)) {
      x->qc[3]++;
      continue;
    }
    if (b->core.qual < 20 || b->core.qual == 255) {
      x->qc[4]++;
      continue;
    }
    const char *name = sam_hdr_tid2name(h, b->core.tid);
    int chrom = 0;
    if (name && strncmp(name, "chr", 3) == 0) {
      char *end;
      chrom = (int)strtol(name + 3, &end, 10);
      if (*end)
        chrom = 0;
    }
    if (chrom < 1 || chrom > 22) {
      x->qc[5]++;
      continue;
    }
    chrom--;
    if (sam_hdr_tid2len(h, b->core.tid) != x->lengths[chrom])
      FAIL("Reference contig length mismatch");
    uint8_t *qs = bam_aux_get(b, "qs");
    double quality = qs ? bam_aux2f(qs) : NAN;
    if (!isfinite(quality) || quality < 10) {
      x->qc[6]++;
      continue;
    }
    if (!bam_aux_get(b, "MM") || !bam_aux_get(b, "ML")) {
      x->qc[7]++;
      continue;
    }
    char *rid = bam_get_qname(b);
    khiter_t k = kh_get(readnames, x->seen, rid);
    if (k != kh_end(x->seen)) {
      x->qc[8]++;
      continue;
    }
    char *saved = strdup(rid);
    if (!saved)
      FAIL("Read ID allocation failed");
    int added;
    k = kh_put(readnames, x->seen, saved, &added);
    if (added < 0) {
      free(saved);
      FAIL("Read ID hash allocation failed");
    }
    x->qc[9]++;
    if (bam_parse_basemod(b, state) < 0)
      FAIL("Malformed MM/ML tags");
    uint32_t *cg = bam_get_cigar(b);
    uint8_t *seq = bam_get_seq(b);
    int q = 0, rev = (flag & BAM_FREVERSE) != 0;
    int64_t rp = b->core.pos;
    uint64_t paired_read = 0;
    for (uint32_t ci = 0; ci < b->core.n_cigar; ci++) {
      int op = bam_cigar_op(cg[ci]), len = bam_cigar_oplen(cg[ci]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        x->qc[10] += len;
        uint32_t j =
            lower(x->sites, x->offsets[chrom], x->offsets[chrom + 1], rp - rev);
        for (;
             j < x->offsets[chrom + 1] && (int64_t)x->sites[j] + rev < rp + len;
             j++) {
          int64_t genomic = (int64_t)x->sites[j] + rev;
          int qp = q + (int)(genomic - rp);
          if (qp < 100 || qp >= b->core.l_qseq - 100 || genomic <= 0 ||
              genomic >= x->lengths[chrom] - 1)
            continue;
          if (bam_seqi(seq, qp) != (rev ? 4 : 2))
            continue; /* G on reverse, C on forward */
          uint32_t *out = x->counts + (size_t)j * 8;
          out[5]++;
          int intact =
              rev ? (bam_seqi(seq, qp - 1) == 2) : (bam_seqi(seq, qp + 1) == 4);
          if (intact)
            out[6]++;
          hts_base_mod mods[16];
          int nm = bam_mods_at_qpos(b, qp, state, mods, 16), qm = -1, qh = -1;
          if (nm < 0 || nm > 16)
            FAIL("Modification parser error/too many modifications");
          for (int z = 0; z < nm; z++)
            if (mods[z].canonical_base == 'C' && mods[z].strand == 0) {
              if (mods[z].modified_base == 'm')
                qm = mods[z].qual;
              if (mods[z].modified_base == 'h')
                qh = mods[z].qual;
            }
          if (qm < 0 || qh < 0) {
            if (intact)
              out[7]++;
            continue;
          }
          if (qm + qh > 256)
            FAIL("Modification probabilities exceed quantization tolerance");
          if (qm + qh > 255)
            x->qc[11]++;
          double pm = (qm + .5) / 256., ph = (qh + .5) / 256.,
                 pc = fmax(0., 1. - pm - ph);
          int win = 0;
          double conf = pc;
          if (pm > conf) {
            win = 1;
            conf = pm;
          }
          if (ph > conf) {
            win = 2;
            conf = ph;
          }
          out[4]++;
          paired_read++;
          if (conf >= .8)
            out[win]++;
          else
            out[3]++;
          uint32_t region = region_offsets[chrom] + x->sites[j] / 100000;
          const double thresholds[4] = {.7, .8, .9, .95};
          for (int ti = 0; ti < 4; ti++)
            if (conf >= thresholds[ti])
              x->regions[(size_t)region * 12 + ti * 3 + win]++;
        }
        q += len;
        rp += len;
      } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP)
        q += len;
      else if (op == BAM_CDEL || op == BAM_CREF_SKIP)
        rp += len;
    }
    if (paired_read)
      x->qc[12]++;
    else
      x->qc[13]++;
  }
  if (ret < -1)
    FAIL("BAM decoding failed before EOF");
  status = 0;
finish:
  if (state)
    hts_base_mod_state_free(state);
  if (b)
    bam_destroy1(b);
  if (h)
    sam_hdr_destroy(h);
  if (f)
    sam_close(f);
  return status;
}
