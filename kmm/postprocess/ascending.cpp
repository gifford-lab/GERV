// Copyright 2014 Daniel Kang

#include <boost/program_options.hpp>

#include <string>

#include "../classes/GradArgs.h"
#include "../util/util.h"

class Ascending: public GradArgs {
 public:
  Ascending() : GradArgs(0, 0) {}

  void load_yall(const char *fname) {
    readd(fname, yall.data(), yall.size());
  }

  void dump_ygrad(const char *fname) {
    ascending_k(KBIG, 1, yall.data(), ygrad.data());
    dumpd(fname, ygrad.data(), ygrad.size());
  }
};

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  po::variables_map vm;

  try {
    desc.add_options()
      ("yout", po::value<std::string>()->required(), "Output directory")
      ("yall", po::value<std::string>()->required(), "Covariates file for heldout.");

    po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
    po::notify(vm);

    fprintf(stdout, "Running with:\nK2: %d\nKBIG: %d\n\n", K2, KBIG);
    boost::program_options::basic_parsed_options<char> parsed_options =
        boost::program_options::parse_command_line(argc, argv, desc);
    for (auto& option : parsed_options.options) {
      std::cout << option.string_key << ": ";
      for (auto& c : option.value)
        std::cout << c;
      std::cout << "\n";
    }
  }
  catch (boost::program_options::required_option& e) {
    fprintf(stdout, "Missing required options\n");
    std::cout << desc << "\n";
    return 1;
  }

  Ascending t;

  t.load_yall(vm["yall"].as<std::string>().c_str());
  t.dump_ygrad(vm["yout"].as<std::string>().c_str());

  return 0;
}
