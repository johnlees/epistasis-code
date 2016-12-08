/*
 * File: common.cpp
 *
 * Program options parsing plus a few other bits
 *
 */

#include "epistasis.hpp"

// Parse command line parameters into usable program parameters
cmdOptions verifyCommandLine(boost::program_options::variables_map& vm, double num_samples)
{
   cmdOptions verified;

   if(vm.count("bacteria"))
   {
      verified.bact_file = vm["bacteria"].as<std::string>();
   }

   if(vm.count("human"))
   {
      verified.human_file = vm["human"].as<std::string>();
   }

   if(vm.count("struct"))
   {
      verified.struct_file = vm["struct"].as<std::string>();
   }

   if(vm.count("chisq"))
   {
      verified.chi_cutoff = stod(vm["chisq"].as<std::string>());
   }

   if (vm.count("pval"))
   {
      verified.log_cutoff = stod(vm["pval"].as<std::string>());
   }

   if (vm.count("chunk_start"))
   {
      verified.chunk_start = vm["chunk_start"].as<long int>();
   }

   if (vm.count("chunk_end"))
   {
      verified.chunk_end = vm["chunk_end"].as<long int>();
   }

   // Error check filtering options
   verified.min_words = 0;
   if (vm.count("mac"))
   {
      int min_ac_in = vm["mac"].as<int>();
      if (min_ac_in >= 0)
      {
         verified.min_ac = min_ac_in;
      }
      else
      {
         throw std::runtime_error("couldn't process min_words argument");
      }
   }
   else
   {
      double maf_in = vm["maf"].as<double>();
      if (maf_in >= 0)
      {
         verified.min_ac = static_cast<unsigned int>(num_samples * maf_in);
      }
      else
      {
         throw std::runtime_error("could process maf argument");
      }
   }

   if (verified.min_ac > num_samples)
   {
      throw std::runtime_error("couldn't process maf/mac argument");
   }
   verified.max_ac = num_samples - verified.min_ac;

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

