#ifndef common_h
#define common_h

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
using namespace arma;

template <typename T>
std::string to_string(T const& value);
void chk_mat(const mat& x, const std::string& varName, const std::string& type);
arma::mat unique_rows(arma::mat x);
double factorial_f(double k);

#endif
