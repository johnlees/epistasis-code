/*
 *
 * seercommon.hpp
 * Header file for seercommon
 * Shared functions between seer and kmds
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <chrono>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <thread>
#include <exception>
#include <sys/stat.h>
#include <regex>

// gzstream headers
#include <gzstream.h>

// Boost headers
#include <boost/program_options.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

// Armadillo/dlib headers
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include <dlib/matrix.h>

// Classes
#include "pair.hpp"

// Constants
const std::string VERSION = "0.1";
//    Default options
const double maf_default = 0.01;
const long int max_length_default = 100;
const double chisq_default = "1";
const double pval_default = "1";

typedef dlib::matrix<double,0,1> column_vector;

// Structs
struct cmdOptions
{
   double log_cutoff;
   double chi_cutoff;

   long int size;
   size_t min_words;
   size_t max_words;

   std::string bact_file;
   std::string human_file;
};

//TODO
// Function headers for each cpp file


