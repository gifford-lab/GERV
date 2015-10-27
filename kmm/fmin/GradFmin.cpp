// Copyright 2014 Daniel Kang

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <thread>
#include <utility>
#include <vector>

#include "../consts.h"
#include "../util/util.h"
#include "./GradFmin.h"

using std::thread;
using std::vector;

GradFmin::GradFmin(GradArgs *m_, GradArgs *h, const int max_it) :
    master(m_), heldout(h), MAX_IT(max_it), iter(0), eta_to_iter() {
}

double GradFmin::minimize() {
  const double start_score = heldout->heldout_feval();
  fprintf(stderr, "Starting score is: %f\n", start_score);
  double start_eta = master->find_eta(KBIG);
  fprintf(stderr, "Starting eta set to: %f\n", start_eta);

  const double mulfact = 0.3;
  double eta = start_eta;
  for (; iter < MAX_IT; ) {
    f_(log(eta));
    eta *= mulfact;
  }

  return 0;
}


double GradFmin::f_(const double x) {
  const double eta = exp(x);

  if (!eta_to_iter.empty()) {
    pair<double, int> hotstart = std::make_pair(FLT_MAX, 0);
    for (auto& p : eta_to_iter)
      if (p.first > eta && p.first < hotstart.first)
        hotstart = p;

    master->load_params(hotstart.second);
  }

  eta_to_iter.push_back(std::make_pair(eta, iter));

  master->update_eta(eta);
  const double prev_heldout = heldout->heldout_feval();
  return master->worker(heldout, iter++, prev_heldout);
}
