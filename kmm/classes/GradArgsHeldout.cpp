// Copyright 2014 Daniel Kang

#include <math.h>

#include <string>
#include <thread>

#include "./GradArgs.h"
#include "../util/util.h"

using std::thread;

double GradArgs::feval_explx0(ftype *explx0, const int64_t s, const int64_t e, const bool pen) {
  int64_t SIZE = e - s;
  const ftype *counts = get_counts(s);
  const gtype *genome = get_genome(s);

  double feval = 0;
  #pragma omp parallel for reduction(+:feval)
  for (int64_t i = 0; i < SIZE; i++) {
    double tmp = 0;
    double rate = explx0[i];
    if (genome[i] < 0)
      tmp = rate - x0[1];
    else
      tmp = rate - x0[0];
    feval += counts[i] * tmp - exp(tmp);
  }

  if (pen) {
    const double eta_eff = eta;
    #pragma omp parallel for reduction(-:feval)
    for (size_t i = 0; i < yall.size(); i++)
      feval -= eta_eff * fabs(yall[i]);
  }

  // We return -feval since we want to cast the problem as minimization.
  return -feval;
}


double GradArgs::heldout_feval(const bool swap_hist, const int64_t EPOCH_SIZE) {
  // double t1 = my_time();
  (void) swap_hist;
  (void) EPOCH_SIZE;

  ftype *explx0 = (ftype*) calloc(ALLOC_SIZE(end - start), sizeof(ftype));
  // For flanks
  explx0 += 2*K;

  sum_explx0(explx0 + K, start + K, end - 2*K);
  double feval = feval_explx0(explx0, start, end - 2*K, false);

  free(explx0 - 2*K);

  // fprintf(stderr, "heldout score: %f, eta: %f, x0: %f, x0nt: %f, time: %f\n",
  //         feval, eta, x0, x0nt, (my_time() - t1)/1000.0);
  return feval;
}


double GradArgs::heldout_feval() {
  return heldout_feval(false, 0);
}
