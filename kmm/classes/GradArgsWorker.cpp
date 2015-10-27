// Copyright 2014 Daniel Kang

#include <float.h>
#include <math.h>

#include <algorithm>
#include <array>
#include <random>
#include <thread>
#include <utility>
#include <vector>

#include "./GradArgs.h"
#include "../util/util.h"

using std::array;
using std::pair;
using std::thread;
using std::vector;

void GradArgs::sum_explx0_grad_worker(ftype *explx0, const int64_t s, const int64_t e) {
  int64_t SIZE = e - s;
  const gtype *genome = get_genome(s);
  for (int64_t i = 0; i < SIZE; i++) {
    int k1 = genome[i];
    if (k1 < 0) continue;

    int ind = k1 * KALLOC;
    for (int64_t j = 0; j< 2*K / RESOL; j++) {
      ftype yg = ygrad[ind + j];
      int64_t yind = i - K + j * RESOL;
      for (int64_t k = 0; k < RESOL; k++)
        explx0[yind + k] += yg;
    }
  }
}
void GradArgs::sum_explx0_grad(ftype *buf, const int64_t s, const int64_t e) {
  const int64_t SUM_SIZE = (e - s) + (3*K);
  const int64_t STEP_SIZE = SUM_SIZE / NUM_THREADS / 2;

  for (int i = 0; i < 2; i++) {
    vector<thread> t(NUM_THREADS);
    for (size_t j = 0; j < t.size(); j++) {
      int64_t s_offset = -K + STEP_SIZE * (2*j + i);
      int64_t e_offset = -K + STEP_SIZE * (2*j + i + 1);
      t[j] = thread(&GradArgs::sum_explx0_grad_worker, this, buf + s_offset,
                    s + s_offset, s + e_offset);
    }
    for (size_t j = 0; j < t.size(); j++)
      t[j].join();
  }
}

void GradArgs::sum_beta(ftype *buf, const int64_t s, const int64_t e) {
  const int64_t SIZE = e - s;
  const ftype *covars = get_covars(s);

  for (int j = 0; j < NUM_COV; j++) {
    #pragma omp parallel for
    for (int64_t i = -K_BETA; i < SIZE + K_BETA; i++)
      for (int ind = 0; ind < K2_BETA; ind++) {
        const int shift = ind - K_BETA;
        buf[i + shift] += beta[j * K2_BETA + ind] * covars[i * NUM_COV + j];
      }
  }
}


void GradArgs::sum_explx0(ftype *buf, const int64_t s, const int64_t e) {
  ascending_k(KBIG, 1, yall.data(), ygrad.data());
  sum_explx0_grad(buf, s, e);
  sum_beta(buf, s, e);
}

void GradArgs::update_x0(ftype *explx0, const int64_t s, const int64_t e) {
  int64_t SIZE = e - s;
  const gtype *genome = get_genome(s);
  const ftype *counts = get_counts(s);

  // Update x0
  double dx0[2], dx0_ = 0, dx0nt = 0;
  double cx0[2], cx0_ = 0, cx0nt = 0;
  #pragma omp parallel for reduction(+:dx0_, dx0nt, cx0_, cx0nt)
  for (int64_t i = 0; i < SIZE; i++) {
    if (genome[i] < 0) {
      dx0_ += exp(explx0[i]);
      cx0_ += counts[i];
    } else {
      dx0nt += exp(explx0[i]);
      cx0nt += counts[i];
    }
  }
  dx0[0] = dx0nt; dx0[1] = dx0_;
  cx0[0] = cx0nt; cx0[1] = cx0_;

  for (size_t i = 0; i < x0.size(); i++) {
    x0[i] = log((dx0[i] + 1) / (cx0[i] + 1));
  }
}

void GradArgs::update_beta(ftype *buf, const ftype eps_beta, const int64_t s, const int64_t e) {
  const int64_t SIZE = e - s;
  const ftype *covars = get_covars(s);

  dbeta.clear();
  #pragma omp parallel for
  for (int ind = 0; ind < K2_BETA; ind++) {
    const int shift = ind - K_BETA;
    for (int j = 0; j < NUM_COV; j++) {
      ftype tmp = 0;
      for (int64_t i = -K_BETA; i <= SIZE + K_BETA; i++)
        tmp += buf[i + shift] * covars[i * NUM_COV + j];
      dbeta[j * K2_BETA + ind] = tmp;
    }
  }


  for (size_t i = 0; i < beta.size(); i++)
    beta[i] += dbeta[i] * eps_beta;
}


void GradArgs::explx0_error(ftype *explx0, const int64_t s, const int64_t e) {
  int64_t SIZE = e - s;
  const gtype *genome = get_genome(s);
  const ftype *counts = get_counts(s);

  #pragma omp parallel for
  for (int64_t i = 0; i < SIZE; i++) {
    ftype rate = explx0[i];
    if (genome[i] < 0)
      explx0[i] = counts[i] - exp(rate - x0[1]);
    else
      explx0[i] = counts[i] - exp(rate - x0[0]);
  }

  for (int64_t i = -K2; i < 0; i++) explx0[i] = 0;
  for (int64_t i = SIZE; i < SIZE + K2; i++) explx0[i] = 0;
}

void GradArgs::compute_gradient(ftype *explx0, const ftype eps_eff,
                                const int64_t s, const int64_t e) {
  ArrayF<YALL_SIZE, 0> yall_ada_tmp;
  ftype **yk_ada_tmp = RENAME(yall_ada_tmp.data(), KBIG, KALLOC);

  int64_t SIZE = e - s;
  const gtype *genome = get_genome(s);

  #pragma omp parallel for
  for (int64_t i = -K; i <= SIZE + K; i++) {
    ftype gtmp = 0;
    for (int64_t j = 0; j < RESOL; j++)
      gtmp += explx0[i + j];
    explx0[i] = gtmp;
  }

  // Update the other parameters.
  ftype* ygrad_tmp = yk_ada_tmp[KBIG];
  #pragma omp parallel for
  for (int64_t i = -K; i <= SIZE + K; i++) {
    int k1 = genome[i];
    if (k1 < 0) continue;

    int ind = k1 * KALLOC;
    for (int64_t j = 0, yind = i - K; j < 2*K/RESOL; j++, yind += RESOL)
#if defined(SAFE_THREADING)
      #pragma omp atomic
#endif
      ygrad_tmp[ind + j] += explx0[yind];
  }
  descending_k(yall_ada_tmp.data(), KBIG, 1, KALLOC);

  const ftype eta_eff = eta;
  #pragma omp parallel for
  for (size_t i = 0; i < yall.size(); i++) {
    const ftype step_size = eps_eff / (kmer_ct[ i / KALLOC ] + 1);
    const ftype primal = yall[i] + yall_ada_tmp[i] * step_size;
    const ftype soft_tr = fabs(primal) - eta_eff * step_size;
    if (soft_tr <= 0)
      yall[i] = 0;
    else
      yall[i] = sign(primal) * soft_tr;
  }

  free(yk_ada_tmp);
}

double GradArgs::get_penalty(ftype *buf) {
  double tmp = 0;
  #pragma omp parallel for reduction(+:tmp)
  for (size_t i = 0; i < YALL_SIZE; i++)
    tmp += fabs(buf[i]);
  return tmp;
}




double GradArgs::find_eps(int64_t s, int64_t e) {
  eps_offset = 0;

  // double eps_eff = 0.000001;
  double eps_eff = 0.0004;

  memset(explx0_main - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));
  memset(explx0_prev - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));

  ArrayF<NUM_BETA, 0> beta_old;
  beta_old.copy(beta);
  xall.copy(yall);

  sum_explx0(explx0_main, s, e);
  double fv0 = feval_explx0(explx0_main, s, e, true);
  fprintf(stderr, "fv0: %f\n", fv0);

  // copy explx0_main to explx0_pgrad since we want to keep explx0_main.
  for (int64_t i = 0; i < ALLOC_SIZE(e - s); i++) {
    explx0_prev[i - K2] = explx0_main[i - K2];
  }

  explx0_error(explx0_prev, s, e);
  update_beta(explx0_prev, eps_eff * beta_mult, s, e);
  compute_gradient(explx0_prev, eps_eff, s, e);

  memset(explx0_prev - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));

  sum_explx0(explx0_prev, s, e);

  double fv_prev = fv0;
  double fv1 = feval_explx0(explx0_prev, s, e, true);
  fprintf(stderr, "fv1: %f\n", fv1);
  while ( fv1 < fv_prev + 1e-3 ) {
    eps_eff *=2;
    #pragma omp parallel for
    for (int64_t i = 0; i < (e - s); i++) {
      double diff = explx0_prev[i] - explx0_main[i];
      explx0_prev[i] += diff;
    }
    fv_prev = fv1;
    fv1 = feval_explx0(explx0_prev, s, e, true);
    fprintf(stderr, "fv1: %f, eps: %f\n", fv1, eps_eff);
  }

  memset(explx0_prev - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));
  memset(explx0_main - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));

  beta.copy(beta_old);
  yall.copy(xall);

  update_eps(eps_eff / 4);
  return fv0;
}


bool GradArgs::params_valid() {
  for (size_t i = 0; i < yall.size(); i++)
    if (!std::isfinite(yall[i]))
      return false;
  for (size_t i = 0; i < x0.size(); i++)
    if (!std::isfinite(x0[i]))
      return false;
  return true;
}


double GradArgs::worker(GradArgs *heldout, const int iter, const double prev_heldout) {
  (void) prev_heldout;

  double prev_fv = find_eps(start + K, end - K) + 0.1;

  vector<double> fpast{prev_fv, prev_fv};

  dump_params(fpast[0], iter);
  int pen_inversions = 0;

  double nest_theta = 1;

  double prev_penalty = 0;
  double penalty = 0;

  // reset explx0
  memset(explx0_prev - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));
  memset(explx0_main - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));

  // init explx0_prev to the previous iter variate:
  sum_explx0(explx0_prev, start + K2, end - K2);

  int t_backtrack = 0;
  double heldout_score = 0;

  beta_nest.copy(beta);
  for (int i = 0; (i < 15 || pen_inversions < 3 || t_backtrack < 5) && i < 50; i++) {
    double t1 = my_time();

    // accumulate current parameters into explx0 and then
    sum_explx0(explx0_main, start + K2, end - K2);

    // eval f(x) at y_k
    double fv = feval_explx0(explx0_main, start + K2, end - K2, true);

    // current variable state:
    // yall = x_k
    // xall = x_{k-1}
    // explx0_main = generated from x_k
    // explx0_prev = generated from x_{k-1}
    // fv = f(x_k)
    // fpast[0] = f(x_{k-1})
    penalty = get_penalty(yall.data());
    prev_penalty = get_penalty(xall.data());

    // if nonmonotone, do backtrcking:
    double alpha = 0.7;
    if ( (fv >= fpast[0] * 1.001 && i < 5) || (fv >= fpast[0] + 1 && i >=5) ) {
      // proposal for eps: backtrack along the arc explx0_main <-> explx0_prev.
      // alpha: backtrack fraction:
      fprintf(stderr, "backtracking: ");
      double eps_propose = eps;
      double fv_back = fv;
      while (fv_back > fpast[0] && eps_propose > 0.00001) {
        // perform backtrack
        eps_propose = eps * alpha;
        #pragma omp parallel for
        for (int64_t j = 0; j < (end - start); j++)
          explx0_main[j] = alpha * explx0_main[j] + (1 - alpha) * explx0_prev[j];
        penalty = alpha * penalty + (1-alpha) * prev_penalty;
        update_x0(explx0_main, start + K2, end - K2);
        fv_back = feval_explx0(explx0_main, start + K2, end - K2, false) + eta * penalty;
        fprintf(stderr, "%f,", fv_back);
      }
      fprintf(stderr, "\n");
      // given eps_propose, reset iterate to x_{k-1} and do the sum_explx0 with eps_propose
      yall.copy(xall);
      beta.copy(beta_nest);
      update_eps(eps_propose / 4);
      nest_theta = 1;
      // since param changed, re-calc explx0
      memset(explx0_main - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));
      sum_explx0(explx0_main, start + K2, end - K2);
      // bookkeeping
      i--;
      t_backtrack = 0;
      pen_inversions = 0;
    } else {  // if monotone update bookkeeing vars (accept fv).
      fpast.insert(fpast.begin(), fv);
      pen_inversions += prev_penalty > penalty;
      t_backtrack++;
      // if we've been successful for a while, propose incrasing stepsize
      // idea: keep increasing stepsize to check what max feasible stepsize would hav been.
      if (t_backtrack > 5) {
        fprintf(stderr, "increasing stepsize: ");
        double eps_propose = eps;
        double fv_prev = fv;
        double fv_new = fv;
        dumpd((out_dir + "explx0_main.bin").c_str(), explx0_main - K2, ALLOC_SIZE(end-start));
        while (fv_new <= fv_prev && eps_propose <= eps * 16) {
          penalty = penalty + (penalty - prev_penalty);
          eps_propose = eps_propose * 2;
          #pragma omp parallel for
          for (int64_t j = 0; j < (end - start); j++) {
            double diff = explx0_main[j] - explx0_prev[j];
            explx0_main[j] = explx0_main[j] + diff;
          }
          fv_prev = fv_new;
          update_x0(explx0_main, start + K2, end - K2);
          fv_new = feval_explx0(explx0_main, start + K2, end - K2, false) + eta * penalty;
          fprintf(stderr, "%f,", fv_new);
        }
        update_eps(std::max(eps_propose / 8.0, eps));
        t_backtrack = 0;
        // since param changed, re-calc explx0
        readd((out_dir + "explx0_main.bin").c_str(), explx0_main - K2, ALLOC_SIZE(end-start));
        // memset(explx0_main - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));
        // sum_explx0(explx0_main, start + K2, end - K2);
      }
      #pragma omp parallel for
      for (int64_t j = 0; j < (end - start); j++)
        explx0_prev[j] = explx0_main[j];
    }

    penalty = get_penalty(yall.data());

    // now perform the update using gradient
    // if reset triggered, we are now at x_k

    // current state of variables
    // yall holds current variable (x_k)
    // xall holds stale variable (x_{k-1})

    double gamma = (nest_theta - 2) / (nest_theta + 1);
    nest_theta++;

    #pragma omp parallel for
    for (size_t j = 0; j < yall.size(); j++) {
      // equivalent to y_k  = x_k + gamma * ( x_k - x_{k-1} )
      double ynext = yall[j] + gamma * (yall[j] - xall[j]);
      // now that y_k is updated, now longer need x_{k-1}, so set xall = x_k
      xall[j] = yall[j];
      // set yall = y_k.
      yall[j] = ynext;
    }
    for (size_t j = 0; j < beta.size(); j++) {
      double bnext = beta[j] + gamma * (beta[j] - beta_nest[j]);
      beta_nest[j] = beta[j];
      beta[j] = bnext;
    }

    update_x0(explx0_main, start + K2, end - K2);

    explx0_error(explx0_main, start + K2, end - K2);
    update_beta(explx0_main, eps * beta_mult, start + K2, end - K2);
    compute_gradient(explx0_main, eps, start + K2, end - K2);

    // clean explx0
    memset(explx0_main - K2, 0, ALLOC_SIZE(end - start) * sizeof(ftype));

    // now
    // yall = x_{k+1}
    // xall = x_k

    heldout_score = heldout->heldout_feval(true, 0);
    fprintf(stderr, "iter: %d fv: %f held: %f x0nt: %f, x0: %f, penalty: %f time: %f\n",
      i, fv, heldout_score, x0[0], x0[1], penalty, (my_time()-t1)/1000.0);
  }

  fprintf(stderr, "finished iteration: %d \n", iter);
  dump_params(heldout_score, iter);
  return fpast[0];
}
