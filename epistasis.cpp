/*
 * File: epistasis.cpp
 *
 * Main control loop for epistasis
 *
 */

#include "epistasis.hpp"

// Constants
const std::string VERSION = "0.1";
//    Default options
const double maf_default = 0.01;
const long int max_length_default = 100;
const std::string chisq_default = "1";
const std::string pval_default = "1";
const double convergence_limit = 10e-8;
const unsigned int max_nr_iterations = 1000;
const double se_limit = 3;

// Starting value for beta vectors (except intercept)
// Should be >0. This value is based on RMS in example study
const double bfgs_start_beta = 1;

int main (int argc, char *argv[])
{
   // Read line of file to get size
   // Read cmd line options and process
   // Open the struct object, check size ok
   // Make a new pair obj to be used throughout, set covars
   // Open the human file, read through until required chunk
   // Read in all the bacterial variants
   // Read a human line (define sub to do this here) - outer loop
   // Set x
   // Set y - inner loop
   // Check maf filter
   // Run chisq test, check filter
   // Run logistic test
   // Print results
   // Set a new y (automatically resets stats)

   // Program description
   std::cerr << "epistasis: pairwise correlations between variants in two populations\n";

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: epistasis -b bacterial_snps_indels.csv.gz -h human_snps.csv.gz --struct all_structure\n\n"
         << "For full option details run seer -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   // Open human file to get number of samples
   std::vector<std::string> human_variant;
   if (vm.count("human"))
   {
      igzstream human_in;
      human_in.open(vm["pheno"].as<std::string>().c_str());
      human_variant = readCsvLine(human_in);
   }
   else
   {
      throw std::runtime_error("--human option is compulsory");
   }

   size_t num_samples = human_variant.size();

   // Error check command line options
   cmdOptions parameters = verifyCommandLine(vm, num_samples);

   // Get mds values
   arma::mat mds;
   int use_mds = 0;
   if (fileStat(parameters.struct_file))
   {
      arma::mat mds;
      mds.load(parameters.struct_file);

      if (mds.n_rows != num_samples)
      {
         throw std::runtime_error("Number of rows in MDS matrix does not match number of samples");
      }
      else
      {
         use_mds = 1;
         std::cerr << "WARNING: Struct file loaded. IT IS UP TO YOU to make sure the order of samples is the same as in each matrix\n";
      }
   }

   // calculate the null log-likelihood
   arma::mat x(num_samples, 1, arma::fill::ones);
   if (use_mds)
   {
      x = join_rows(x, mds);
   }

   // Open the human variant ifstream, and read through until the required
   // block is reached
   igzstream human_file;
   human_file.open(parameters.human_file.c_str());

   if (parameters.chunk_start != 0 && parameters.chunk_end != 0)
   {
      if (parameters.chunk_start >= parameters.chunk_end)
      {
         throw std::runtime_error("chunk start greater than or equal to chunk end");
      }
      else
      {
         std::cerr << "Reading to chunk position: line " << parameters.chunk_start << std::endl;
         std::string line;
         for (int i = 0; i < parameters.chunk_start - 1; i++)
         {
            std::getline(human_file, line);
         }
      }
   }

   // Read in all the bacterial variants (3Mb compressed - shouldn't be too bad
   // in this form I hope)
   std::cerr << "Reading in all bacterial variants" << std::endl;
   igzstream bacterial_file;
   bacterial_file.open(parameters.bact_file.c_str());

   std::vector<Pair> all_pairs;
   long int line_nr = 1;
   while (bacterial_file)
   {
      std::vector<std::string> bacterial_variant = readCsvLine(bacterial_file);
      Pair bact_in(num_samples);
      bact_in.add_y(bacterial_variant, line_nr);
      if (use_mds)
      {
         bact_in.add_covar(mds);
      }

      all_pairs.push_back(bact_in);
      line_nr++;
   }

   // Write a header
   std::cerr << "Starting association tests" << std::endl;
   std::string header = "human_line\tbact_line\thuman_af\tbacterial_af\tchisq_p_val\tlogistic_p_val\tlrt_p_val\tbeta\tcomments";
   std::cout << header << std::endl;

   // Read the block of human lines
   human_file.open(parameters.human_file.c_str());

   line_nr = 1;
   long int read_pairs = 0;
   long int tested_pairs = 0;
   long int significant_pairs = 0;
   while (human_file)
   {
      std::vector<std::string> human_variant = readCsvLine(human_file);
      // Test each human variant against every bacterial variant
      for (auto it = all_pairs.begin(); it < all_pairs.end(); it++)
      {
         it->add_x(human_variant, line_nr);

         // maf filter
         std::tuple<double,double> mafs = it->maf();
         if (!(std::get<0>(mafs) < parameters.min_af || std::get<0>(mafs) > parameters.max_af || std::get<1>(mafs) < parameters.min_af || std::get<1>(mafs) > parameters.max_af))
         {
            it->chisq_p(chiTest(*it));
            read_pairs++;

            if (it->chisq_p() < parameters.chi_cutoff)
            {
               doLogit(*it);

               // Likelihood ratio test
               double null_ll = nullLogLikelihood(it->get_x(), it->get_y());
               it->p_val(likelihoodRatioTest(*it, null_ll));

               tested_pairs++;
               if (it->p_val() < parameters.log_cutoff)
               {
                  significant_pairs++;
               }
            }

            std::cout << *it << std::endl;
         }
      }

      if (line_nr > (parameters.chunk_end - parameters.chunk_start))
      {
         break;
      }
      else
      {
         line_nr++;
      }
   }

   std::cerr << "Processed " << line_nr * all_pairs.size() << " total pairs. Of these:\n";
   std::cerr << "\tPassed maf filter:\t" << read_pairs << std::endl;
   std::cerr << "\tPassed chi^2 filter:\t" << tested_pairs << std::endl;
   std::cerr << "\tPassed p-val (logistic) filter:\t" << significant_pairs << std::endl;
   std::cerr << "Done.\n";
}

std::vector<std::string> readCsvLine(std::istream& is)
{
   std::vector<std::string> variant;
   std::string line;
   std::getline(is, line);

   std::stringstream line_stream(line);
   std::string value;

   while(std::getline(line_stream, value, ','))
   {
      variant.push_back(value);
   }
   return variant;
}

