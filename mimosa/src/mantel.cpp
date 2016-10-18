#include <Rcpp.h>
// #include <algorithm>    // std::shuffle
// #include <random>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {

  // clone a into b to leave a alone
  Rcpp::NumericVector b = Rcpp::clone(a);

  std::random_shuffle(b.begin(), b.end(), randWrapper);

  return b;
}

//' Make matrix of permutations for mantel test with a row for each permutation and a column for each sample
//'
//' @param nSamples number of samples
//' @param nPerm number of permutations
// [[Rcpp::export]]
NumericMatrix make_perm_mat(double nSamples, double nPerm) {
  NumericMatrix out(nPerm, nSamples);
  NumericVector SampSeq(nSamples);
  for(int j = 0; j < nSamples; j++){
    SampSeq[j] = j+1;
  }
  for (int i = 0; i < nPerm; i++){
    NumericMatrix::Row zzrow = out(i, _);
    zzrow = randomShuffle(SampSeq);
  }
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// /*** R
// timesTwo(42)
// randomShuffle(1:14)
// make_perm_mat(4, 6)
// */
