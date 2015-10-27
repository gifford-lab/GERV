// Copyright 2014 Daniel Kang

#include <utility>
#include <vector>

#include "../consts.h"
#include "../classes/GradArgs.h"
#include "../util/util.h"

#ifndef FMIN_GRADFMIN_H_
#define FMIN_GRADFMIN_H_

using std::pair;
using std::vector;

typedef double (*func_eval)(double, void *);

class GradFmin {
  GradArgs *master;
  GradArgs *heldout;
  const int MAX_IT;
  int iter;
  vector<pair<double, int> > eta_to_iter;

 public:
  GradFmin(GradArgs *m_, GradArgs *h_, const int max_it);

  ~GradFmin() {}
  GradFmin(const GradFmin&) = delete;
  GradFmin& operator=(const GradFmin&) = delete;

  void common_run();

  // Function evaluation
  double f_(const double x);

  double minimize();
};

#endif  // FMIN_GRADFMIN_H_
