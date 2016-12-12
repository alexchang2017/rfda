#include "common.h"

//' trapz: trapezoidal rule to approximate the integral values
//'
//' Returns approximation of integral.
//'
//' @param m A numeric matrix to be divided into list of matrices.
//' @param margin The margin of the matrix to split.
//' @param f A integer vector to split the matrix.
//' @return A list of matrices.
//' @examples
//' x <- matrix(rnorm(30), 6, 8)
//' splitMat(x, 1, rep(1:3, each = 2))
//' splitMat(x, 2, rep(1:4, each = 2))
//' @export
// [[Rcpp::export]]
Rcpp::List splitMat(const arma::mat& m, const double& margin, const arma::vec& f) {
  chk_mat(m, "m", "double");
  chk_mat(f, "f", "double");
  if (std::abs(margin - floor(margin)) > 1e-6 || (margin != 1.0 && margin != 2.0))
    Rcpp::stop("margin must be 1 or 2.");
  if (margin == 1.0 && m.n_rows != f.n_elem) {
    Rcpp::stop("The length of f must be equal to the number of rows of m when margin = 1.");
  } else if (margin == 2.0 && m.n_cols != f.n_elem) {
    Rcpp::stop("The length of f must be equal to the number of columns of m when margin = 2.");
  }

  vec uni_f = unique(f);
  mat tmp;
  Rcpp::List outList(uni_f.n_elem);
  if (margin == 1) {
    for (uword i = 0; i < uni_f.n_elem; i++)
    {
      tmp = m.rows(find(f == uni_f(i)));
      outList[i] = Rcpp::wrap(tmp);
    }
  } else if (margin == 2) {
    for (uword i = 0; i < uni_f.n_elem; i++)
    {
      tmp = m.cols(find(f == uni_f(i)));
      outList[i] = Rcpp::wrap(tmp);
    }
  }
  return outList;
}
