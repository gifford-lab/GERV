// Copyright 2014 Daniel Kang

#include "./Loader.h"

inline int get_kmer(const char *genome, const int64_t ind, const int ksize) {
  int kmer = 0;
  for (int i = 0; i < ksize; i++)
    kmer |= ((int)genome[ind+i]) << (2*i);
  return kmer;
}


inline bool ind_ok(const char *genome, const int64_t ind, const int ksize) {
  bool ok = 1;
  for (int j = 0; j < ksize; j++)
    ok &= (genome[ind+j] < 4 && genome[ind+j] >= 0);
  return ok;
}


gtype* Loader::load_genome(const int64_t SIZE, const int KS,
                           const int64_t offset, const char *fname) {
  gtype *genome = (gtype *) calloc(ALLOC_SIZE(SIZE), sizeof(gtype));

  char *buf = (char *) fread_c(offset, sizeof(char), SIZE, fname);

  #pragma omp parallel for
  for (int64_t i = 0; i < SIZE; i++) {
    if (ind_ok(buf, i, KS))
      genome[i] = get_kmer(buf, i, KS);
    else
      genome[i] = -1;
  }

  free(buf);

  return genome;
}


ftype* Loader::load_counts(const int64_t SIZE, const int max_reads, const int window_size,
                           const int64_t offset, const char *fname) {
  ftype* counts = (ftype*) fread_c(offset, sizeof(ftype), SIZE, fname);

  // Truncate the read counts
  int64_t num_truncated_bases = 0;
  int64_t total_truncated_reads = 0;
  #pragma omp parallel for reduction(+:num_truncated_bases, total_truncated_reads)
  for (int64_t i = 0; i < SIZE; i++) {
    if (counts[i] > max_reads) {
      num_truncated_bases++;
      total_truncated_reads += counts[i] - max_reads;
      counts[i] = max_reads;
    }
  }

  fprintf(stderr, "num bases: %ld bases cut: %ld reads cut: %ld\n",
          SIZE, num_truncated_bases, total_truncated_reads);

  // Smooth the counts
#if defined(SAFE_THREADING)
  #pragma omp parallel for
#endif
  for (int64_t i = 0; i < SIZE - window_size; i++) {
    double tmp = 0;
    for (int j = 0; j < window_size; j++)
      tmp += counts[i + j];
    counts[i] = tmp / window_size;
  }

  return counts;
}
