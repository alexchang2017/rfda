#include "common.h"
#include "locPoly1d.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

// a parallel worker to compute quantiles of y give bandwidth along xout
struct Worker_quantData: public RcppParallel::Worker {
  const double& bandwidth;
  const vec& probs;
  const vec& x;
  const vec& y;
  const vec& w;
  const vec& xout;
  vec& new_x;
  mat& new_y;
  vec& new_w;

  Worker_quantData(const double& bandwidth, const vec& probs, const vec& x, const vec& y, const vec& w,
                   const vec& xout, vec& new_x, mat& new_y, vec& new_w):
    bandwidth(bandwidth), probs(probs), x(x), y(y), w(w), xout(xout), new_x(new_x), new_y(new_y), new_w(new_w) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (uword i = begin; i < end; ++i)
    {
      // find the data between xout[i] - bandwidth and xout[i] + bandwidth
      uvec in_windows = linspace<uvec>(0, x.n_elem-1, x.n_elem);
      in_windows = in_windows.elem(find(all(join_rows(x <= xout(i) + bandwidth,
                                                      x >= xout(i) - bandwidth), 1)));
      if (in_windows.n_elem > 0)
      {
        // find the median of x in the windows
        new_x(i) = median(x.elem(in_windows));
        // find the median of w in the windows
        new_w(i) = median(w.elem(in_windows));
        // find the quantiles of y in the windows
        mat ly = y.rows(in_windows);
        new_y.row(i) = quantileCpp(ly, probs).t();
      } else
      {
        // return nan if there is no data in the windows
        new_x(i) = datum::nan;
        new_w(i) = datum::nan;
        new_y.row(i).fill(datum::nan);
      }
    }
  }
};

//' One-dimensional kernel local polynominal smoother of quantiles
//'
//' Perform one-dimensional kernel local polynominal smoother of quantiles corresponding to the
//' given probabilities for data \code{(x,y)} with weight \code{w} on \code{xout}.
//'
//' @param bandwidth A single numerical value. The kernel smoothing parameter.
//' @param probs A numeric vector with values between 0 and 1. The probabilities of quantiles.
//' @param x A vector, the variable of of x-axis.
//' @param y A vector, the variable of of y-axis. \code{y[i]} is corresponding value of \code{x[i]}.
//' @param w A vector, the weight of data. \code{w[i]} is corresponding value of \code{x[i]}.
//' @param xout A vector, vector of output time points. It should be a sorted vecotr.
//' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
//' @param drv An integer, the order of derivative.
//' @param degree An integer, the degree of polynomial.
//' @return A estimated value on \code{xout} by one-dimensional kernel local polynominal smoother of quantiles
//' corresponding to the given probabilities.
//' @references
//' \enumerate{
//'   \item https://www.r-statistics.com/2010/04/quantile-loess-combining-a-moving-quantile-window-with-loess-r-function
//' }
//' @examples
//' N <- 100
//' x <- runif(N, 0, 10)
//' y <- rnorm(N)
//' xout <- sort(runif(200, 0, 10))
//' est <- locPoly1d(1.2, x, y, rep(1, N), xout, 'gauss', 0, 1)
//' require(pipeR)
//' require(data.table)
//' est_quant <- locQuantPoly1d(1.2, c(0.25, 0.5, 0.75), x, y, rep(1, N), xout, 'gauss', 0, 1) %>>%
//'   data.table %>>% setnames(c("Q25", "Q50", "Q75"))
//' linesDF <- est_quant[ , `:=`(x = xout, loess = est)] %>>% melt.data.table("x", value.name = "y")
//' require(ggplot2)
//' ggplot(data.frame(x,y), aes(x, y)) + geom_point() +
//'   geom_line(aes(x, y, colour = variable), data = linesDF) +
//'   labs(colour = "Type")
//' @export
// [[Rcpp::export]]
arma::mat locQuantPoly1d(const double& bandwidth, const arma::vec& probs, const arma::vec& x,
                         const arma::vec& y, const arma::vec& w, const arma::vec& xout,
                         const std::string& kernel, const double& drv, const double& degree){
  // check data
  chk_mat(probs, "probs", "double");
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  chk_mat(w, "w", "double");
  if (x.n_elem != y.n_elem || x.n_elem != w.n_elem)
    Rcpp::stop("The lengths of x, y and w must be equal.\n");
  chk_mat(xout, "xout", "double");
  if (!xout.is_sorted())
    Rcpp::stop("xout must be sorted.\n");
  if (!is_finite(bandwidth) || bandwidth <= 0.0)
    Rcpp::stop("bandwidth must be a positive number.\n");
  if (any(probs > 1.0) || any(probs < 0.0))
    Rcpp::stop("probs must be a numeric vector with values between 0 and 1.\n");

  // remove the data with 0 weights
  uvec actObs = find(w != 0);
  vec xw = x.elem(actObs), yw = y.elem(actObs), ww = w.elem(actObs);
  // allocate output data
  mat new_y = zeros<mat>(xout.n_elem, probs.n_elem);
  vec new_x = zeros<vec>(xout.n_elem), new_w = zeros<vec>(xout.n_elem);

  // compute parallely the quantiles given bandwidth
  Worker_quantData quantData(bandwidth, probs, xw, yw, ww, xout, new_x, new_y, new_w);
  RcppParallel::parallelFor(0, xout.n_elem, quantData);
  // remove the infinite values
  uvec idxNonfinte = find_nonfinite(new_w);
  if (idxNonfinte.n_elem > 0){
    if ((double) idxNonfinte.n_elem / (double) new_w.n_elem > 0.25)
      Rcpp::stop("The bandwidth is too small, please pick another one!");
    new_w.elem(idxNonfinte).zeros();
    new_x.elem(idxNonfinte).zeros();
    new_y.rows(idxNonfinte).zeros();
  }

  // smooth the quantiles with local kernel polynominal smoother
  mat est = zeros<mat>(xout.n_elem, probs.n_elem);
  for (uword i = 0; i < probs.n_elem; i++){
    vec tmp_y = new_y.col(i);
    est.col(i) = locPoly1d_cpp(bandwidth, new_x, tmp_y, new_w, xout, kernel, drv, degree);
  }

  return est;
}
