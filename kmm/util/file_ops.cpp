// Copyright 2014 Daniel Kang

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "./util.h"

void* fread_c(size_t start, size_t size, size_t num_el, const char *fname) {
  void *tmp = calloc(ALLOC_SIZE(num_el), size);
  FILE *fin = fopen(fname, "r");
  fseek(fin, start, SEEK_SET);
  size_t read_size = fread(tmp, size, num_el, fin);
  fclose(fin);

  if (read_size != num_el) {
    fprintf(stderr, "failed to read file %s\n", fname);
    exit(1);
  }

  return tmp;
}


void dumpd(const char *fname, const ftype *data, size_t n) {
  FILE *fout = fopen(fname, "w");
  fwrite(data, sizeof(ftype), n, fout);
  fclose(fout);
}

void readd(const char *fname, ftype *data, size_t n) {
  FILE *fin = fopen(fname, "r");
  size_t size = fread(data, sizeof(ftype), n, fin);
  fclose(fin);
  if (size != n) {
    fprintf(stderr, "failed to read file");
    exit(1);
  }
}


void dump_out(const char *out_dir, const char *var_name, const ftype *out, int size, int iter) {
  char out_fname[1024];

  snprintf(out_fname, sizeof(out_fname), "%s/%s_%d.bin", out_dir, var_name, iter);
  dumpd(out_fname, out, size);
}


void load_out(const char *out_dir, const char *var_name, ftype *out, int size, int iter) {
  char out_fname[1024];

  snprintf(out_fname, sizeof(out_fname), "%s/%s_%d.bin", out_dir, var_name, iter);
  readd(out_fname, out, size);
}

