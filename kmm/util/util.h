// Copyright 2013 Daniel Kang

#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

#include "../consts.h"

#ifndef UTIL_UTIL_H_
#define UTIL_UTIL_H_

inline int get_num_cores() {
  return sysconf(_SC_NPROCESSORS_ONLN);
}

void error_and_quit(const char *error);

void check_mem(void *arr, const char *arr_name);

void check_thread(int rc, char *err_msg);

void* fread_c(size_t start, size_t size, size_t num_el, const char *fname);

void dumpd(const char *fname, const ftype *data, size_t n);

void readd(const char *fname, ftype *data, size_t n);

void dump_out(const char *out_dir, const char *var_name, const ftype *out, int size, int iter);

void load_out(const char *out_dir, const char *var_name, ftype *out, int size, int iter);

inline ftype sign(ftype x) {
  return (0 < x) - (x < 0);
}

inline ftype thresh(ftype x, ftype eps, ftype eta) {
  if (fabs(x) < eps*eta)
    return 0;
  else
    return sign(x)*(fabs(x) - eps*eta);
}

void threshx(ftype *out, const ftype *in, const ftype eps,
             const ftype eta, const int size);

ftype** RENAME(const ftype *input, const int MAXK, const int stride);
void ascending_k(int MAXK, int MINK, const ftype *karr_flat, ftype *out);
void descending_k(ftype *karr, int MAXK, int MINK, int stride);


inline void copy_prev(ftype *prev, const ftype *cur, const int64_t size) {
  memcpy(prev, cur, sizeof(ftype)*size);
}

inline bool file_exists(const char *fname) {
  struct stat buffer;
  return stat(fname, &buffer) == 0;
}

inline bool filesize_valid(const char *fname, int64_t min_filesize) {
  FILE *f = fopen(fname, "r");
  fseek(f, 0L, SEEK_END);
  // If the stream can read past SEEK_END it usually means the file is special.
  if (fgetc(f) != EOF)
    return true;
  return ftell(f) >= min_filesize;
}

double my_time();

#endif  // UTIL_UTIL_H_
