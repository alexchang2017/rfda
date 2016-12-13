#include "common.h"
#include <cmath>

// function to print message on R console
void RMessage(const std::string& msg){
  Rcpp::Function RMessage("message");
  RMessage(msg);
}

// function to check whether the input data with correct type
void chk_mat(const mat& x, const std::string& varName, const std::string& type){
  if (!is_finite(x))
    Rcpp::stop(varName + " must be numerical.\n");
}

// function to compute factorial in C++ implementation
// [[Rcpp::export]]
double factorial_f(double k){
  if (std::abs(k - floor(k)) > 1e-6)
    Rcpp::stop("factorial_f does not support float input.");

  if (k == 0) {
    return 1.0;
  } else if (k > 0) {
    double out = 1;
    for (double i = k; i > 0; --i)
      out *= i;
    return out;
  } else {
    return datum::nan;
  }
}

// function to compute kernel density with 4 kinds of kernel
arma::vec kernelDensity(const arma::vec& x, const std::string& kernel){
  vec kx(x.n_elem);
  if (kernel == "epan") {
    kx = (1 - square(x)) * 0.75;
    kx.elem(find(abs(x) >= 1)).zeros();
  } else if (kernel == "gaussvar") {
    kx = exp(-0.5*square(x)) / sqrt(2 * datum::pi) % (1.25 - 0.25 * square(x));
  } else if (kernel == "quar") {
    kx = square(1 - square(x)) * (15.0 / 16.0);
    kx.elem(find(abs(x) >= 1)).zeros();
  }  else if (kernel == "gauss") {
    kx = exp(-0.5*square(x)) / sqrt(2 * datum::pi);
  }
  return kx;
}

// function to compute quantiles with interpolation method in C++ implementation
// [[Rcpp::export]]
arma::vec quantileCpp(const arma::vec& x, const arma::vec& probs){
  // check data
  chk_mat(x, "x", "double");
  chk_mat(probs, "probs", "double");
  if (any(probs < 0) || any(probs > 1))
    Rcpp::stop("There is values in probs outside [0, 1].");

  // sort data
  vec x_sort = sort(x);
  // find the location of quantiles
  vec index = (x.n_elem - 1) * probs;
  // find the floor and ceiling integer index
  uvec lo = conv_to<uvec>::from(floor(index)), hi = conv_to<uvec>::from(ceil(index));
  // find the weight to interpolate
  vec h = (index - conv_to<vec>::from(lo));
  // interpolate the quantiles
  vec qs = (1 - h) % x_sort.elem(lo) + h % x_sort.elem(hi);
  return qs;
}

// function to compute trapezoidal numerical integration
// [[Rcpp::export]]
arma::mat trapz_cpp(const arma::vec& x, const arma::mat& y){
  // check data
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  if (y.n_rows != x.n_elem)
    Rcpp::stop("The number of rows or length of y must be equal to the length of x.");
  if (x.n_elem <= 1)
    Rcpp::stop("The length of x must be greater than 2.");

  return diff(x).t() * (y.rows(0, y.n_rows-2) + y.rows(1, y.n_rows-1)) * 0.5;
}
