/*
 * pair.hpp
 * Header file for Pair class
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <cstdlib>
#include <string>
#include <iterator>
#include <vector>
#include <exception>

// Armadillo/dlib headers
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

std::string pair_comment_default = "NA";

class Pair
{
   public:
      // Initialisation
      Pair(int number_samples);

      // nonmodifying operations
      long int bact_line() const { return _bact_line; }
      long int human_line() const { return _human_line; }

      std::tuple<double,double> maf() const { return std::make_tuple (_maf_x, _maf_y); }
      double chisq_p() const { return _chisq_p; }
      double p_val() const { return _lrt_p; }
      double beta() const { return _beta; }
      double se() const { return _se; }
      std::string comments() const { return _comment; }
      int firth() const { return _firth; }

      arma::vec get_x() const { return _x; }
      arma::vec get_y() const { return _y; }
      arma::mat get_x_design(); const { return join_rows(_x, _covars); }

      // Modifying operations
      void set_line_nr(const long int line_nr) { _line_nr = line_nr; }
      void p_val(const double pvalue) { _lrt_p = pvalue; }
      void chisq_p(const double pvalue) { _chisq_p = pvalue; }
      void beta(const double b) { _beta = b; }
      void standard_error(const double se) { _se = se; }
      void set_maf(const double maf) { _maf = maf; }
      void firth(const int set_firth) { _firth = set_firth; }

      void add_comment(const std::string& new_comment); // this is defined in pair.cpp
      void add_x(const std::vector<std::string>& variant); // this is defined in pair.cpp
      void add_y(const std::vector<std::string>& variant); // this is defined in pair.cpp
      void add_covar(const arma::mat& covars); // this is defined in pair.cpp
      void reset_stats(); // this is defined in pair.cpp

   private:
      int _number_samples;

      long int _bact_line;
      long int _human_line;

      arma::vec _y;
      arma::mat _x;
      arma::mat _covars;
      int _covars_set;

      double _maf_x;
      double _maf_y;
      double _chisq_p;
      double _lrt_p;
      double _log_likelihood;
      double _beta;
      double _se;
      std::string _comment;

      int _firth;
};

// Overload output operator
std::ostream& operator<<(std::ostream &os, const Significant_kmer& k);

