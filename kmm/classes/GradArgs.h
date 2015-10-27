// Copyright 2014 Daniel Kang

#include <boost/program_options.hpp>

#include <string>

#include "../consts.h"
#include "./ArrayF.h"
#include "./Loader.h"

#ifndef CLASSES_GRADARGS_H_
#define CLASSES_GRADARGS_H_

using boost::program_options::variables_map;


class GradArgs : public Loader {
  std::string out_dir = "/mnt/";
  int64_t eps_offset = -1;
  const double beta_mult;

 protected:
  static ArrayF<YALL_SIZE, 0> yall;
  static ArrayF<YALL_SIZE, 0> xall;
  static ArrayF<NUM_BETA, 0> beta;
  static ArrayF<NUM_BETA, 0> beta_nest;
  static ArrayF<NUM_BETA, 0> dbeta;
  static ArrayF<YALL_SIZE / KALLOC, 0> kmer_ct;
  static ftype **kmer_ctk;
  static ftype **yk;

  static ftype *explx0_main;
  static ftype *explx0_prev;

  static ArrayF<YGRAD_SIZE, 0> ygrad;

  static ArrayF<2, 0> x0;
  static ArrayF<2, 0> x0_prev;

  static ftype epoch;
  static double eta;
  static double eps;

  static double fv_epoch;

 public:
  // Worker constructor
  GradArgs(const int64_t N, const variables_map& vm);
  // Heldout constructor.
  GradArgs(const variables_map& vm) : Loader(vm), beta_mult(0) {}
  // Don't do anything
  GradArgs(const int64_t s, const int64_t e) : Loader(s, e), beta_mult(0) {}
  // Only load genome, covars
  GradArgs(const int64_t N, const int64_t offset,
           const std::string& genome_fname, const std::string& covars_fname) :
      Loader(N, offset, genome_fname, covars_fname), beta_mult(0) {}

  virtual ~GradArgs() {}

  GradArgs(const GradArgs&) = delete;
  GradArgs& operator=(const GradArgs&) = delete;

  void optimize_beta(ftype *explx0, const int64_t s, const int64_t e);
  void scatter_data();

  double find_eta(const int ksize);

  void update_eta(double e);
  void update_eps(double e);

  void dump_params(const ftype fval, const int iter);
  void load_params(const int iter);

  double get_penalty(ftype *yall);

  double feval_explx0(ftype *buf, const int64_t s, const int64_t e, const bool pen);

  void sum_explx0_grad_worker(ftype *explx0, const int64_t s, const int64_t e);
  void sum_explx0_grad(ftype *explx0, const int64_t s, const int64_t e);
  void sum_explx0(ftype *buf, const int64_t s, const int64_t e);
  void sum_beta(ftype *buf, const int64_t s, const int64_t e);

  void update_x0(ftype *buf, const int64_t s, const int64_t e);
  void update_beta(ftype *buf, const ftype eps_beta, const int64_t s, const int64_t e);

  void explx0_error(ftype *buf, const int64_t s, const int64_t e);

  void set_eps_offset();
  double find_eps(int64_t s, int64_t e);

  bool params_valid();

  void compute_gradient(ftype *buf, const ftype eps_eff, const int64_t s, const int64_t e);
  double worker(GradArgs *heldout, const int iter, const double prev_heldout);

  double heldout_feval(const bool swap_hist, const int64_t EPOCH_SIZE);
  double heldout_feval();
};

#endif  // CLASSES_GRADARGS_H_
