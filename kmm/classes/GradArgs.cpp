// Copyright 2014 Daniel Kang

#include <math.h>

#include <gsl/gsl_multimin.h>

#include <algorithm>
#include <string>
#include <vector>

#include "./GradArgs.h"
#include "../util/util.h"

using std::vector;


ArrayF<YALL_SIZE, 0> GradArgs::yall;
ArrayF<YALL_SIZE, 0> GradArgs::xall;
ArrayF<NUM_BETA, 0> GradArgs::beta;
ArrayF<NUM_BETA, 0> GradArgs::beta_nest;
ArrayF<NUM_BETA, 0> GradArgs::dbeta;
ArrayF<YALL_SIZE / KALLOC, 0> GradArgs::kmer_ct;
ftype** GradArgs::kmer_ctk = RENAME(kmer_ct.data(), KBIG, 1);
ftype** GradArgs::yk = RENAME(yall.data(), KBIG, KALLOC);

ftype* GradArgs::explx0_main;
ftype* GradArgs::explx0_prev;

ArrayF<YGRAD_SIZE, 0> GradArgs::ygrad;

ArrayF<2, 0> GradArgs::x0;
ArrayF<2, 0> GradArgs::x0_prev;

ftype GradArgs::epoch = 1;
double GradArgs::eta = 0;
double GradArgs::eps = 0;

double GradArgs::fv_epoch = 0;


GradArgs::GradArgs(const int64_t N, const variables_map& vm) :
    Loader(N, vm), beta_mult(1.0 / N) {
  double t1 = my_time();

  out_dir = vm["out_dir"].as<std::string>();
  epoch = 1;

  explx0_main = (ftype*) calloc(ALLOC_SIZE(end - start), sizeof(ftype));
  explx0_main += K2;

  explx0_prev = (ftype*) calloc(ALLOC_SIZE(end - start), sizeof(ftype));
  explx0_prev += K2;

  optimize_beta(explx0_main, start + K2, end - K2);

  fprintf(stderr, "Total setup time: %f\n", (my_time() - t1)/1000.0);
}


void GradArgs::optimize_beta(ftype *explx0, const int64_t s, const int64_t e) {
  const int64_t SIZE = e - s;
  const double eps_beta = beta_mult * 15;  // TODO(ddkang): Fix eps
  const double tts = my_time();
  fprintf(stderr, "eps_beta: %f\n", eps_beta);

  const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
  gsl_multimin_fdfminimizer *solver;
  gsl_vector *x = gsl_vector_alloc(NUM_BETA);
  for (size_t i = 0; i < NUM_BETA; i++) gsl_vector_set(x, i, 0.0);
  gsl_multimin_function_fdf my_func;

  struct SolverData {
    friend class GradArgs;
    GradArgs *ga;
    ftype *explx0_;
    int64_t SIZE_, s_, e_;
    void clear_explx0() {
      memset(explx0_ - K2, 0, (SIZE_ + K2) * sizeof(ftype));
    }
  } SD = {this, explx0, SIZE, s, e};

  my_func.n = NUM_BETA;
  my_func.f = [](const gsl_vector *v, void *params) {
    SolverData *data = (SolverData*) params;
    for (size_t i = 0; i < NUM_BETA; i++) data->ga->beta[i] = gsl_vector_get(v, i);
    data->clear_explx0();
    data->ga->sum_beta(data->explx0_, data->s_, data->e_);
    data->ga->update_x0(data->explx0_, data->s_, data->e_);
    return data->ga->feval_explx0(data->explx0_, data->s_, data->e_, false);
  };
  my_func.df = [](const gsl_vector *v, void *params, gsl_vector *df) {
    SolverData *data = (SolverData*) params;
    for (size_t i = 0; i < NUM_BETA; i++) data->ga->beta[i] = gsl_vector_get(v, i);
    data->clear_explx0();
    data->ga->sum_beta(data->explx0_, data->s_, data->e_);
    data->ga->update_x0(data->explx0_, data->s_, data->e_);
    data->ga->explx0_error(data->explx0_, data->s_, data->e_);
    data->ga->update_beta(data->explx0_, 0, data->s_, data->e_);
    for (size_t i = 0; i < NUM_BETA; i++) gsl_vector_set(df, i, -data->ga->dbeta[i]);
  };
  my_func.fdf = [](const gsl_vector *v, void *params, double *f, gsl_vector *df) {
    SolverData *data = (SolverData*) params;
    for (size_t i = 0; i < NUM_BETA; i++) data->ga->beta[i] = gsl_vector_get(v, i);
    data->clear_explx0();
    data->ga->sum_beta(data->explx0_, data->s_, data->e_);
    data->ga->update_x0(data->explx0_, data->s_, data->e_);
    const double fv = data->ga->feval_explx0(data->explx0_, data->s_, data->e_, false);
    data->ga->explx0_error(data->explx0_, data->s_, data->e_);
    data->ga->update_beta(data->explx0_, 0, data->s_, data->e_);
    *f = fv;
    for (size_t i = 0; i < NUM_BETA; i++) gsl_vector_set(df, i, -data->ga->dbeta[i]);
  };
  my_func.params = &SD;

  solver = gsl_multimin_fdfminimizer_alloc(T, NUM_BETA);

  gsl_multimin_fdfminimizer_set(solver, &my_func, x, eps_beta * 100, 100);
  for (int iter = 0; iter < 500; iter++) {
    const double t1 = my_time();

    int status = gsl_multimin_fdfminimizer_iterate(solver);
    if (status) break;

    status = gsl_multimin_test_gradient(solver->gradient, 1e-3);
    if (status == GSL_SUCCESS) fprintf(stderr, "Minimum found at:\n");
    fprintf(stderr, "beta iter: %d, fv: %f, time: %f\n", iter, solver->f, (my_time()-t1)/1000.0);
  }

  for (size_t i = 0; i < NUM_BETA; i++) beta[i] = gsl_vector_get(solver->x, i);

  fprintf(stderr, "total beta time: %f\n", (my_time()-tts)/1000.0);
}


double GradArgs::find_eta(const int ksize) {
  assert(ksize <= KBIG);
  DEFINE_CONSTS(ksize);

  const gtype *genome = get_genome(0);

  vector<ftype> ckmer(XSIZE_ALL / KALLOC);
  ftype **ckmer_k = RENAME(&ckmer[0], ksize, 1);


  const int mask = (1 << (2*ksize)) - 1;
  // Generate kmer counts
  for (int64_t i = start; i < end; i++) {
    int k1 = genome[i];
    if (k1 < 0) continue;
    kmer_ctk[ksize][k1 & mask]++;
  }

  // Compute the expected counts.
  memset(explx0_main - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));
  sum_explx0(explx0_main, start + K2, end - K2);

  // Find the correct eta
  const double ex0 = exp(-x0[0]);
  for (int64_t i = start; i < end; i++) {
    if (genome[i] < 0) continue;
    for (int j = 0; j < RESOL; j++)
      ckmer_k[ksize][genome[i] & mask] += exp(explx0_main[i + j] - x0[0]);
  }

  double max_eta = 0;
  for (int kmer = 0; kmer < (1 << (2*ksize)); kmer++) {
    double kmer_count = kmer_ctk[ksize][kmer];
    double kmer_reads = ckmer_k[ksize][kmer];
    double eta_tmp = fabs(kmer_reads - kmer_count * ex0);
    max_eta = std::max(eta_tmp, max_eta);
  }

  descending_k(kmer_ct.data(), KBIG, 1, 1);

  return max_eta;
}


void GradArgs::update_eta(double e) {
  eta = e;
}


void GradArgs::update_eps(double e) {
  fprintf(stderr, "eps set to: %f\n", e);
  eps = e;
}


// THIS DUMPS HISTORY
void GradArgs::dump_params(const ftype fval, const int iter) {
  dump_out(out_dir.c_str(), "xall", xall.data(), xall.size(), iter);
  dump_out(out_dir.c_str(), "yall", yall.data(), yall.size(), iter);
  dump_out(out_dir.c_str(), "beta", beta.data(), beta.size(), iter);
  dump_out(out_dir.c_str(), "beta_nest", beta_nest.data(), beta_nest.size(), iter);
  dump_out(out_dir.c_str(), "x0", x0.data(), x0.size(), iter);
  dump_out(out_dir.c_str(), "x0_prev", x0_prev.data(), x0_prev.size(), iter);
  const ftype eta_tmp = eta;
  dump_out(out_dir.c_str(), "eta", &eta_tmp, 1, iter);
  dump_out(out_dir.c_str(), "heldout", &fval, 1, iter);
}


void GradArgs::load_params(const int iter) {
  load_out(out_dir.c_str(), "xall", xall.data(), xall.size(), iter);
  load_out(out_dir.c_str(), "yall", yall.data(), yall.size(), iter);
  load_out(out_dir.c_str(), "beta", beta.data(), beta.size(), iter);
  load_out(out_dir.c_str(), "beta_nest", beta_nest.data(), beta_nest.size(), iter);
  load_out(out_dir.c_str(), "x0", x0.data(), x0.size(), iter);
  load_out(out_dir.c_str(), "x0_prev", x0_prev.data(), x0_prev.size(), iter);
}
