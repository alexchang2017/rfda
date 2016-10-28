#ifndef common_h
#define common_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
using namespace arma;

void RMessage(const std::string& msg);
template <typename T>
std::string to_string(T const& value);
void chk_mat(const mat& x, const std::string& varName, const std::string& type);
arma::mat unique_rows(arma::mat x);
double factorial_f(double k);
arma::vec kernelDensity(const arma::vec& xin, const std::string& kernel);
arma::vec quantileCpp(const arma::vec& x, const arma::vec& probs);
arma::mat trapz_cpp(const arma::vec& x, const arma::mat& y);
#endif
