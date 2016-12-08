/*
 * File: stats
 *
 * Implements chi^2 test
 *
 */

#include "epistasis.hpp"

const double normalArea = pow(2*M_PI, -0.5);

// Basic chi^2 test, using contingency table
double chiTest(Pair& p)
{
   arma::mat x = p.get_x();
   arma::vec y = p.get_y();

   double chisq = 0;

   // Contigency table
   //          human 0   human 1   human 2
   // bact 0   a         b         c
   // bact 1   d         e         f
   //
   // Use doubles for compatibility with det function in arma::mat
   double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;

   arma::vec::const_iterator j = y.begin();
   for (arma::vec::const_iterator i = x.begin(); i!=x.end(); ++i)
   {
      if (*j == 0) {
         if (*i == 0){
            a++;
         } else if (*i == 1){
            b++;
         } else {
            c++;
         }
      } else {
         if (*i == 0){
            d++;
         } else if (*i == 1){
            e++;
         } else {
            f++;
         }
      j++;
   }

   arma::mat::fixed<2, 3> table = {a, b, c, d, e, f};
#ifdef SEER_DEBUG
   arma::Mat<int>::fixed<2, 3> tab_out = {int (a), int (b), int (c), int(d), int(e), int(f)};
   std::cerr << tab_out << "\n";
#endif

   int N = accu(table);

   if (N == 0)
   {
      throw std::logic_error("Empty table for chisq test\n");
   }

   // Treat as invalid if any entry is 0 or 1, or if more than one entry < 5
   // Mark as needing to use Firth regression
   int low_obs = 0;
   for (auto obs = table.begin(); obs != table.end(); ++obs)
   {
      if (*obs <= 1 || (*obs <= 5 && ++low_obs > 2))
      {
         p.add_comment("bad-chisq");
         p.firth(1);
         break;
      }
   }

   // With Yates' continuity correction
   double chisq = 0;

   double total = accu(table);
   arma::vec col_sum = sum(table, 0);
   arma::vec row_sum = sum(table, 1);

   for (int i = 0; i < 2; ++i)
   {
      for (int j = 0; j < 3; j++)
      {
         double expected = row_sum[i] * col_sum[j] / total;
         chisq += pow((std::abs(table[i,j] - expected) - 0.5), 2) / expected;
      }
   }

   boost::math::chi_squared_distribution chi_dist(2);
   double p_value = 1 - boost::math::cdf(chi_dist, chisq);

   if (p_value = 0)
   {
      p_value = normalPval(pow(chisq, 0.5));
      p.add_comment("chi-large");
   }
#ifdef SEER_DEBUG
   std::cerr << "chisq:" << chisq << "\n";
   std::cerr << "chisq p: " << p_value << "\n";
#endif
   return p_value;
}

// Fit null models for null log-likelihoods
double nullLogLikelihood(const arma::mat& x, const arma::vec& y)
{
   double null_ll = 0;

   if (x.n_cols > 1)
   {
      Pair null_pair(x.n_rows);

      null_pair.add_x(x, 0);
      null_pair.add_y(y, 0);

      doLogit(null_pair);
      null_ll = null_pair.log_likelihood();
   }
   else
   {
      dlib::matrix<double,1,1> intercept;

      intercept(0) = log(mean(y)/(1-mean(y))); // null is: intercept = log-odds of success
      LogitLikelihood likelihood_fit(x, y);
      null_ll = likelihood_fit(intercept);
   }
   return null_ll;
}

// Likelihood-ratio test
double likelihoodRatioTest(Pair& p, const double null_ll)
{
   double log_likelihood = p.log_likelihood();
   double lrt_p = 1;
   if (log_likelihood == 0 || null_ll == 0)
   {
      k.add_comment("zero-ll");
   }
   else
   {
      lrt = pow(2*(log_likelihood - null_ll), 0.5);

      if (lrt > 0)
      {
         lrt_p = normalPval(lrt);
      }
   }
   return lrt_p;
}

// Returns p-value for a test statistic that is >0 and standard normally distributed
double normalPval(double testStatistic)
{
   double p_val = 0;
   if (testStatistic < 5)
   {
      boost::math::normal s;

      p_val = 2 * (1 - boost::math::cdf(s, testStatistic));
   }
   else
   {
      // For large z need to use a bound
      // See http://stats.stackexchange.com/questions/13690/how-to-compute-the-probability-associated-with-absurdly-large-z-scores
      //
      // Upper bound
      // S(z) <= phi(z)/z
      // cdf = 1-(0.5 * S(z))
      // At z = 5 correct to +/- 2.5%
#ifdef SEER_DEBUG
      std::cerr << "using erfc bound rather than 'exact' function\n";
#endif
      p_val = 2 * exp(-0.5*pow(testStatistic,2))*normalArea/testStatistic;
   }

   return p_val;
}


