#include "trapz.h"

// [[Rcpp::export]]
arma::mat trapz_cpp(const arma::vec& x, const arma::mat& y){
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  if (y.n_rows != x.n_elem)
    Rcpp::stop("The number of rows or length of y must be equal to the length of x.");
  if (x.n_elem <= 1)
    Rcpp::stop("The length of x must be greater than 2.");

  return diff(x).t() * (y.rows(0, y.n_rows-2) + y.rows(1, y.n_rows-1)) * 0.5;
}
