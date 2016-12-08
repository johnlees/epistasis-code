/*
 * File: common.cpp
 *
 * Program options parsing plus a few other bits
 *
 */

#include "epistasis.hpp"

// Parse command line parameters into usable program parameters
cmdOptions verifyCommandLine(boost::program_options::variables_map& vm, const std::vector<Sample>& samples)
{
   cmdOptions verified;

   verified.max_length = vm["max_length"].as<long int>();

   if(vm.count("kmers"))
   {
      verified.kmers = vm["kmers"].as<std::string>();
   }

   if(vm.count("output"))
   {
      verified.output = vm["output"].as<std::string>();
   }

   if(vm.count("chisq"))
   {
      verified.chi_cutoff = stod(vm["chisq"].as<std::string>());
   }

   if (vm.count("pval"))
   {
      verified.log_cutoff = stod(vm["pval"].as<std::string>());
   }

   // Verify MDS options in a separate function
   // This is pc, size and number of threads
   verifyMDSOptions(verified, vm);

   verified.filter = 1;
   if (vm.count("no_filtering"))
   {
      verified.filter = 0;
   }
   else
   {
      // Error check filtering options
      verified.min_words = 0;
      if (vm.count("min_words"))
      {
         int min_words_in = vm["min_words"].as<int>();
         if (min_words_in >= 0)
         {
            verified.min_words = min_words_in;
         }
         else
         {
            badCommand("min_words", std::to_string(min_words_in));
         }
      }
      else
      {
         double maf_in = vm["maf"].as<double>();
         if (maf_in >= 0)
         {
            verified.min_words = static_cast<unsigned int>(samples.size() * maf_in);
         }
         else
         {
            badCommand("maf", std::to_string(maf_in));
         }
      }

      if (verified.min_words > samples.size())
      {
         badCommand("min_words/maf", std::to_string(verified.min_words));
      }
      verified.max_words = samples.size() - verified.min_words;
   }

   verified.print_samples = 0;
   if (vm.count("print_samples"))
   {
      verified.print_samples = 1;
   }

   return verified;
}

// Conversion functions required as code is a mix of dlib and armadillo
// matrices
// This could obviously be improved...
arma::vec dlib_to_arma(const column_vector& dlib_vec)
{
   arma::vec converted(dlib_vec.nr());

   for (unsigned int i = 0; i < dlib_vec.nr(); ++i)
   {
      converted(i) = dlib_vec(i);
   }

   return converted;
}

column_vector arma_to_dlib(const arma::vec& arma_vec)
{
   column_vector converted;
   converted.set_size(arma_vec.n_elem);

   for (unsigned int i = 0; i < arma_vec.n_elem; ++i)
   {
      converted(i) = arma_vec(i);
   }

   return converted;
}

// Inverts a symmetric positive matrix, checking for errors
// Not passed by ref, creates a copy. Right thing to do?
arma::mat inv_covar(arma::mat A)
{
   // Try the default. Internally this uses Cholesky decomposition and back
   // solves. For large condition numbers it fails.
   arma::mat B;
   if (!inv_sympd(B, A))
   {
      // If the Cholesky decomposition fails, try pseudo-inverse
      // This uses SVD:
      // A = U*S*V.t() => A^-1 = V*S^-1*U.t()
      // and ignores small values in the S matrix
      if (!arma::pinv(B, A))
      {
         std::cerr << "A matrix inversion failed!" << std::endl;
      }
   }

   return B;
}

// Check for file existence
int fileStat(const std::string& filename)
{
   struct stat buffer;
   int success = 1;

   if (stat (filename.c_str(), &buffer) != 0)
   {
      std::cerr << "Can't stat input file: " << filename << "\n";

      success = 0;
   }

   return success;
}

