#include "common.h"
#include <cmath>
#include <sstream>

template <typename T>
std::string to_string(T const& value) {
  std::stringstream sstr;
  sstr << value;
  return sstr.str();
}

// the function to check the input data type
void chk_mat(const mat& x, const std::string& varName, const std::string& type){
  if (!is_finite(x))
    Rcpp::stop(varName + " must be a number.\n");

  if (type == "integer") {
    if (all(all(abs(x - floor(x)) < 1e-6)))
      Rcpp::stop(varName + " must be a nonnegative number.\n");
  }
}

// [[Rcpp::export]]
double factorial_f(double k){
  if (std::abs(k - floor(k)) > 1e-6)
    Rcpp::stop("factorial_f does not support float input.");

  if (k == 0) {
    return 1.0;
  } else if (k > 0) {
    double out = 1;
    for (double i = k; i > 0; i--)
      out *= i;
    return out;
  } else {
    return datum::nan;
  }
}
