// Copyright 2014 Daniel Kang

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <boost/program_options.hpp>

#include <string>

#include "../consts.h"

int64_t N, NUM_CHR;

void proc_chr(ftype *out, FILE *fin, int64_t offset, int64_t max_ct) {
  char buf[5000];
  int64_t loc;

  char *ret_val;
  while (!feof(fin)) {
    memset(buf, 0, sizeof(buf));
    ret_val = fgets(buf, sizeof(buf), fin);
    if (ret_val == NULL || strlen(buf) < 2) break;

    sscanf(buf, "%ld", &loc);

    loc += offset - 1;
    if (loc > N) {
      printf("uhoh!");
      exit(1);
    }
    if (out[loc] < max_ct)
      out[loc]++;
  }
}

int main(int argc, char **argv) {
  FILE *fin, *fout, *foffsets;
  const char *in_dir;
  int64_t max_ct;

  namespace po = boost::program_options;
  po::options_description desc("Options");
  po::variables_map vm;

  desc.add_options()
    ("out_file", po::value<std::string>()->required(), "Output file")
    ("reads_dir", po::value<std::string>()->required(), "Directory of reads")
    ("offsets_file", po::value<std::string>()->required(), "File with offsets, one per chromosome.")
    ("max_ct", po::value<int64_t>()->required(), "Maximum count per base")
    ("num_bases", po::value<int64_t>()->required(), "Number of bases")
    ("num_chr", po::value<int64_t>()->required(), "Number of chromosomes");

  try {
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    fout = fopen(vm["out_file"].as<std::string>().c_str(), "w");
    in_dir = vm["reads_dir"].as<std::string>().c_str();
    foffsets = fopen(vm["offsets_file"].as<std::string>().c_str(), "r");
    max_ct = vm["max_ct"].as<int64_t>();
    N = vm["num_bases"].as<int64_t>();
    NUM_CHR = vm["num_chr"].as<int64_t>();
  }
  catch (boost::program_options::required_option& e) {
    fprintf(stdout, "Missing required options\n");
    std::cout << desc;
    return 1;
  }

  fprintf(stderr, "Number of bases: %ld\n", N);
  ftype *out = (ftype*) calloc(N + 10000, sizeof(ftype));
  if (out == NULL) {
    fprintf(stderr, "Failed to allocate out\n");
    exit(1);
  }

  for (int i = 1; i <= NUM_CHR; i++) {
    char buf[5000];
    int64_t offset;
    if (fscanf(foffsets, "%ld\n", &offset) != 1) {
      fprintf(stderr, "Offsets not formatted correctly at %d\n", i);
      exit(1);
    }
    snprintf(buf, sizeof(buf), "%s/allreads-%d.csv", in_dir, i);
    fin = fopen(buf, "r");

    proc_chr(out, fin, offset, max_ct);

    fclose(fin);
    fprintf(stdout, "Finished reading chromosome %d\n", i);
  }

  fwrite(out, sizeof(ftype), N, fout);

  fclose(fout);

  return 0;
}
