// Copyright 2013 Daniel Kang

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include <boost/program_options.hpp>

#include <string>

#include "./consts.h"
#include "classes/GradArgs.h"
#include "util/util.h"

#ifndef MPI_MOTIF_H_
#define MPI_MOTIF_H_

using boost::program_options::variables_map;


void check_po_validity(const variables_map& vm) {
  if (!file_exists(vm["genome"].as<std::string>().c_str()))
    error_and_quit("genome file does not exist.\n");
  if (!file_exists(vm["reads"].as<std::string>().c_str()))
    error_and_quit("reads file does not exist.\n");
  if (!file_exists(vm["covariates"].as<std::string>().c_str()))
    error_and_quit("covariates file does not exist.\n");
  if (!file_exists(vm["heldout_reads"].as<std::string>().c_str()))
    error_and_quit("heldout_reads file does not exist.\n");
  if (!file_exists(vm["heldout_covariates"].as<std::string>().c_str()))
    error_and_quit("heldout_covariates file does not exist.\n");

  if (vm["num_bases"].as<int64_t>() <= 0)
    error_and_quit("num_bases cannot be <= 0\n");
  if (vm["num_bases"].as<int64_t>() > vm["heldout_start"].as<int64_t>())
    error_and_quit("heldout_start must be >= num_bases\n");
  if (vm["heldout_start"].as<int64_t>() < 0)
    error_and_quit("heldout_start cannot be < 0\n");
  if (vm["heldout_size"].as<int64_t>() < 0)
    error_and_quit("heldout_size cannot be < 0\n");
  if (vm["read_max"].as<int64_t>() <= 0)
    error_and_quit("read_max must be > 0\n");
  if (vm["smooth_window_size"].as<int64_t>() <= 0)
    error_and_quit("smooth_window_size must be > 0\n");

  if (!filesize_valid(vm["genome"].as<std::string>().c_str(),
                      vm["num_bases"].as<int64_t>()))
    error_and_quit("File genome must be at least num_bases long.\n");
  if (!filesize_valid(vm["reads"].as<std::string>().c_str(),
                      vm["num_bases"].as<int64_t>() * sizeof(ftype)))
    error_and_quit("File reads must be at least num_bases*sizeof(ftype) long.\n");
  if (!filesize_valid(vm["covariates"].as<std::string>().c_str(),
                      vm["num_bases"].as<int64_t>() * sizeof(ftype) * NUM_COV))
    error_and_quit("File covariates must be at least NUM_COV*num_bases*sizeof(ftype) long.\n");
  if (!filesize_valid(vm["heldout_reads"].as<std::string>().c_str(),
                      vm["heldout_size"].as<int64_t>() * sizeof(ftype)))
    error_and_quit("File heldout_reads is too small.\n");
  if (!filesize_valid(vm["heldout_covariates"].as<std::string>().c_str(),
                      vm["heldout_size"].as<int64_t>() * sizeof(ftype) * NUM_COV))
    error_and_quit("File heldout_covariates is too small.\n");
}

#endif  // MPI_MOTIF_H_
