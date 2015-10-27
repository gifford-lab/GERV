// Copyright 2013 Daniel Kang

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include <string>
#include <vector>

#include "../classes/GradArgs.h"
#include "../util/util.h"

using std::vector;

class Validation : public GradArgs {
 public:
  Validation(const int64_t N, const int64_t offset,
             const std::string& genome_fname, const std::string& covars_fname) :
      GradArgs(N, offset, genome_fname, covars_fname) {}

  void set_ygrad(ftype *tmp, size_t size) {
    memcpy(ygrad.data(), tmp, size * sizeof(ftype));
  }

  void set_beta(ftype *tmp) {
    memcpy(beta.data(), tmp, NUM_BETA * sizeof(ftype));
  }

  void set_x0(ftype x0_, ftype x0nt) {
    x0[1] = x0_;
    x0[0] = x0nt;
  }

  void expxc_worker(const char *out_name) {
    const int64_t SIZE = end - start;
    const gtype *genome = get_genome(0);
    vector<ftype> explx0_v(ALLOC_SIZE(SIZE));
    ftype *explx0 = &explx0_v[0] + 2*K2;

    sum_explx0_grad(explx0 + K, start + K, end - 2*K);
    sum_beta(explx0 + K, start + K, end - 2*K);

    #pragma omp parallel for
    for (int64_t i = 0; i < SIZE; i++) {
      double rate = explx0[i];
      if (genome[i] < 0)
        explx0[i] = exp(rate - x0[1]);
      else
        explx0[i] = exp(rate - x0[0]);
    }

    dumpd(out_name, explx0, SIZE);
  }
};


int main(int argc, char **argv) {
  DEFINE_CONSTS(KBIG);
  int64_t start, end;

  namespace po = boost::program_options;
  po::options_description desc("Options");
  po::variables_map vm;

  desc.add_options()
    ("output", po::value<std::string>()->required(), "Output file")
    ("genome", po::value<std::string>()->required(), "Genome file")
    ("covariates", po::value<std::string>()->required(), "Covariates file.")
    ("xopt", po::value<std::string>()->required(), "xopt file")
    ("beta", po::value<std::string>()->required(), "beta parameters")
    ("x0", po::value<double>()->required(), "x0")
    ("x0nt", po::value<double>()->required(), "x0nt")
    ("start", po::value<int64_t>()->required(), "Starting point in genome")
    ("end", po::value<int64_t>()->required(), "Ending point in genome");

  try {
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    start = vm["start"].as<int64_t>();
    end = vm["end"].as<int64_t>();
    assert(end > start);
    assert(start >= 0);
  }
  catch (boost::program_options::required_option& e) {
    fprintf(stdout, "Missing required options\n");
    std::cout << desc;
    return 1;
  }

  // Get arguments
  ftype* XOPT_ALL = (ftype*) calloc(XSIZE_ALL, sizeof(ftype));
  ftype* XOPT = (ftype*) calloc(XSIZE_GRAD, sizeof(ftype));
  ftype* BETA = (ftype*) calloc(NUM_BETA, sizeof(ftype));

  check_mem(XOPT_ALL, "XOPT");

  readd(vm["xopt"].as<std::string>().c_str(), XOPT_ALL, XSIZE_ALL);
  readd(vm["beta"].as<std::string>().c_str(), BETA, NUM_BETA);
  ascending_k(KBIG, 1, XOPT_ALL, XOPT);

  // Get the estimated rates.
  Validation gargs(end - start, start,
                   vm["genome"].as<std::string>(), vm["covariates"].as<std::string>());

  gargs.set_ygrad(XOPT, XSIZE_GRAD);
  gargs.set_beta(BETA);
  gargs.set_x0(vm["x0"].as<double>(), vm["x0nt"].as<double>());
  gargs.expxc_worker(vm["output"].as<std::string>().c_str());
}
