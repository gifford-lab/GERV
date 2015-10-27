// Copyright 2014 Daniel Kang

#include <omp.h>
#include <boost/program_options.hpp>

#include <string>

#include "../consts.h"
#include "../util/util.h"

using boost::program_options::variables_map;

#ifndef CLASSES_LOADER_H_
#define CLASSES_LOADER_H_

class Loader {
 protected:
  const int NUM_THREADS;
  const int64_t start, end;

 private:
  const gtype *genome_;
  const ftype *counts_;
  const ftype *covars_;

 protected:
  // Don't do anything
  Loader(const int64_t s, const int64_t e) :
      NUM_THREADS(omp_get_max_threads()),
      start(s), end(e), genome_(NULL), counts_(NULL), covars_(NULL) {}

  // Only load genome, covars
  Loader(const int64_t N, const int64_t offset,
         const std::string& genome_fname, const std::string& covars_fname) :
      NUM_THREADS(omp_get_max_threads()),
      start(0), end(N), genome_(load_genome(end - start, KBIG, offset, genome_fname.c_str())),
      counts_(NULL),
      covars_(load_counts((end - start) * NUM_COV, 1, 1, offset * NUM_COV * sizeof(ftype), covars_fname.c_str())) {}

  // Normal constructor
  Loader(const int64_t N, const variables_map& vm) :
      NUM_THREADS(omp_get_max_threads()),
      start(0), end(N),
      genome_(load_genome(end - start, KBIG, 0, vm["genome"].as<std::string>().c_str())),
      counts_(load_counts(end - start, vm["read_max"].as<int64_t>(),
                         vm["smooth_window_size"].as<int64_t>(), 0,
                         vm["reads"].as<std::string>().c_str())),
      covars_(load_counts((end - start) * NUM_COV, 1, 1, 0,
                          vm["covariates"].as<std::string>().c_str())) {}

  // Heldout constructor
  explicit Loader(const variables_map& vm) :
      NUM_THREADS(omp_get_max_threads()),
      start(0), end(vm["heldout_size"].as<int64_t>()),
      genome_(load_genome(end - start, KBIG, vm["heldout_start"].as<int64_t>(),
                         vm["genome"].as<std::string>().c_str())),
      counts_(load_counts(end - start, vm["read_max"].as<int64_t>(),
                         vm["smooth_window_size"].as<int64_t>(), 0,
                         vm["heldout_reads"].as<std::string>().c_str())),
      covars_(load_counts((end - start) * NUM_COV, 1, 1, 0,
                          vm["heldout_covariates"].as<std::string>().c_str())) {}

 public:
  static gtype* load_genome(const int64_t SIZE, const int KS, const int64_t offset,
                            const char *fname);

  static ftype* load_counts(const int64_t SIZE, const int max_reads, const int window_size,
                            const int64_t offset, const char *fname);

  const gtype* get_genome(int64_t i) { return genome_ + i; }
  const ftype* get_counts(int64_t i) { return counts_ + i; }
  const ftype* get_covars(int64_t i) { return covars_ + i; }

  virtual ~Loader() {
    if (genome_ != NULL) free((gtype*) genome_);
    if (counts_ != NULL) free((ftype*) counts_);
  }

  Loader(const Loader&) = delete;
  Loader& operator=(const Loader&) = delete;
};

#endif  // CLASSES_LOADER_H_
