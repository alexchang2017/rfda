#include "common.h"
#include <cmath>
#include <sstream>

void RMessage(const std::string& msg){
  Rcpp::Function RMessage("message");
  RMessage(msg);
}

template <typename T>
std::string to_string(T const& value) {
  std::stringstream sstr;
  sstr << value;
  return sstr.str();
}

// the function to check the input data type
void chk_mat(const mat& x, const std::string& varName, const std::string& type){
  if (!is_finite(x))
    Rcpp::stop(varName + " must be numerical.\n");

  if (type == "integer") {
    if (all(all(abs(x - floor(x)) < 1e-6)))
      Rcpp::stop(varName + " must be integer.\n");
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

arma::vec kernelDensity(const arma::vec& xin, const std::string& kernel){
  vec kx(xin.n_elem);
  if (kernel == "epan")
  {
    kx = (1 - square(xin)) * 0.75;
    kx.elem(find(abs(xin) >= 1)).zeros();
  } else if (kernel == "gaussvar")
  {
    kx = exp(-0.5*square(xin)) / sqrt(2 * datum::pi) % (1.25 - 0.25 * square(xin));
  } else if (kernel == "quar")
  {
    kx = square(1 - square(xin)) * (15.0 / 16.0);
    kx.elem(find(abs(xin) >= 1)).zeros();
  }  else if (kernel == "gauss")
  {
    kx = exp(-0.5*square(xin)) / sqrt(2 * datum::pi);
  }
  return kx;
}
