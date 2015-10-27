// Copyright 2013 Daniel Kang

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "./util.h"

double my_time() {
  /* Linux */
  struct timeval tv;
  gettimeofday(&tv, NULL);

  uint64_t ret = tv.tv_usec;

  /* Adds the seconds (10^0) after converting them to microseconds (10^-6) */
  ret += (tv.tv_sec * 1000000);

  return ret / 1000.0;
}


void error_and_quit(const char *error) {
  fprintf(stderr, "%s\n", error);
  exit(EXIT_FAILURE);
}

void check_mem(void *arr, const char *arr_name) {
  if (arr == NULL) {
    fprintf(stderr, "can't allocate %s\n", arr_name);
    exit(1);
  }
}

void check_thread(int rc, char *err_msg) {
  if (rc != 0) {
    fprintf(stderr, "%s\n", err_msg);
    exit(1);
  }
}


void threshx(ftype *out, const ftype *in, const ftype eps, const ftype eta, const int size) {
  int i;
  for (i = 0; i < size; i++)
    out[i] = thresh(in[i], eps, eta);
}


ftype** RENAME(const ftype *input, const int MAXK, const int stride) {
  ftype** ret = (ftype**) calloc(MAXK + 1, sizeof(ftype*));
  ret[0] = const_cast<ftype*>(input);
  int pow4 = 1;
  for (int k = 1; k <= MAXK; k++) {
    ret[k] = ret[k-1] + (pow4 * stride);
    pow4 *= 4;
  }
  return ret;
}


// Generates dy
void descending_k(ftype *karr_flat, int MAXK, int MINK, int stride) {
  ftype **karr = RENAME(karr_flat, MAXK, stride);

  for (int k = MAXK - 1; k >= MINK; k--) {
    #pragma omp parallel for if (k > 3)
    for (int kmer = 0; kmer < (1 << (k*2)); kmer++) {
      for (int ind = 0; ind < stride; ind++) {
        karr[k][kmer * stride + ind] =
          karr[k+1][(kmer + (0 << (2*k))) * stride + ind] +
          karr[k+1][(kmer + (1 << (2*k))) * stride + ind] +
          karr[k+1][(kmer + (2 << (2*k))) * stride + ind] +
          karr[k+1][(kmer + (3 << (2*k))) * stride + ind];
      }
    }
  }

  free(karr);
}

// Generates y
void ascending_k(int MAXK, int MINK, const ftype *karr_flat, ftype *out) {
  ftype **karr = RENAME(karr_flat, MAXK, KALLOC);
  const int M_MAXK = 1 << (MAXK*2);
  const int out_size = M_MAXK * KALLOC;

  memset(out, 0, out_size * sizeof(ftype));

  for (int k = MINK; k <= MAXK; k++) {
    const int mask = (1 << (2*k)) - 1;
    #pragma omp parallel for
    for (int kmer = 0; kmer < M_MAXK; kmer++)
      for (int ind = 0; ind < KALLOC; ind++)
        out[kmer * KALLOC + ind] += karr[k][(kmer & mask) * KALLOC + ind];
  }

  free(karr);
}

