// Copyright 2013 Daniel Kang

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <boost/program_options.hpp>

#include <string>

#include "./mpi_motif.h"
#include "fmin/GradFmin.h"

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  po::variables_map vm;

  try {
    desc.add_options()
      ("out_dir", po::value<std::string>()->required(), "Output directory")
      ("genome", po::value<std::string>()->required(), "Genome file")
      ("reads", po::value<std::string>()->required(),
       "File containing the numbers of reads at each base")
      ("covariates", po::value<std::string>()->required(),
       "Flat file of covariates. NOTE THE INPUT FORMAT VERY CAREFULLY.")
      ("num_bases", po::value<int64_t>()->required(), "Number of bases.")
      ("heldout_start", po::value<int64_t>()->required(), "Heldout starting position.")
      ("heldout_size", po::value<int64_t>()->required(), "Heldout size.")
      ("read_max", po::value<int64_t>()->required(), "Max reads per base.")
      ("cov_max", po::value<int64_t>()->required(), "THIS OPTION IS IGNORED.")
      ("smooth_window_size", po::value<int64_t>()->required(), "Smoothing window size.")
      ("heldout_reads", po::value<std::string>()->required(), "Reads file for heldout.")
      ("heldout_covariates", po::value<std::string>()->required(), "Covariates file for heldout.");

    po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
    po::notify(vm);

    check_po_validity(vm);
    fprintf(stderr, "Running with:\n");
    fprintf(stderr, "RESOL: %d\n", RESOL);
    fprintf(stderr, "K2: %d\n", K2);
    fprintf(stderr, "KBIG: %d\n", KBIG);
    fprintf(stderr, "NUM_THREADS: %d\n", omp_get_max_threads());
    boost::program_options::basic_parsed_options<char> parsed_options =
        boost::program_options::parse_command_line(argc, argv, desc);
    for (auto& option : parsed_options.options) {
      std::cerr << option.string_key << ": ";
      for (auto& c : option.value)
        std::cerr << c;
      std::cerr << "\n";
    }
  }
  catch (boost::program_options::required_option& e) {
    fprintf(stderr, "Missing required options\n");
    std::cerr << desc << "\n";
    return 1;
  }

  GradArgs *master, *heldout;

  double t1 = my_time();
  master = new GradArgs(vm["num_bases"].as<int64_t>(), vm);
  heldout = new GradArgs(vm);
  fprintf(stderr, "setup time: %f\n", (my_time() - t1)/1000.0);

  GradFmin *minimizer;
  minimizer = new GradFmin(master, heldout, NUM_IT);

  minimizer->minimize();

  return 0;
}
