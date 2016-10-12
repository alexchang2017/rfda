#include "common.h"
#include "interp.h"
#include <sstream>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct Worker_locLinear2d_gauss: public RcppParallel::Worker {
  const vec& bandwidth;
  const mat& xw;
  const vec& yw;
  const vec& ww;
  const vec& countw;
  const vec& out1;
  const vec& out2;
  const std::string& kernel;
  mat& est;

  Worker_locLinear2d_gauss(const vec& bandwidth, const mat& xw, const vec& yw, const vec& ww, const vec& countw,
                           const vec& out1, const vec& out2, const std::string& kernel, mat& est):
    bandwidth(bandwidth), xw(xw), yw(yw), ww(ww), countw(countw), out1(out1), out2(out2),
    kernel(kernel), est(est) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (uword i = begin; i < end; ++i)
    {
      // create model matrix
      mat dx = join_rows(ones<vec>(xw.n_rows), xw);
      // compute demean column
      dx.col(2) -= out2(i);
      for (uword j = 0; j < out1.n_elem; ++j)
      {
        // compute demean column
        dx.col(1) -= out1(j);
        // compute weights
        vec w = ww % exp(-0.5 * square(dx.col(1) / bandwidth(0))) %
          exp(-0.5 * square(dx.col(2) / bandwidth(1))) / (2 * datum::pi);
        if (kernel == "gaussvar")
          w %= (1.25 - 0.25 * square(dx.col(1) / bandwidth(0))) %
            (1.5 - 0.5 * square(dx.col(2) / bandwidth(1)));

        mat dxw = dx;
        dxw.each_col() %= (w % countw);
        vec p = pinv(dxw.t() * dx) * dx.t() * (w % yw);
        est(i, j) = p(0);
        dx.col(1) += out1(j);
      }
    }
  }
};

struct Worker_locLinear2d_nongauss: public RcppParallel::Worker {
  const vec& bandwidth;
  const mat& xw;
  const vec& yw;
  const vec& ww;
  const vec& countw;
  const vec& out1;
  const vec& out2;
  const std::string& kernel;
  mat& est;
  umat& flag;

  Worker_locLinear2d_nongauss(const vec& bandwidth, const mat& xw, const vec& yw, const vec& ww,
                              const vec& countw, const vec& out1, const vec& out2, const std::string& kernel,
                              mat& est, umat& flag):
    bandwidth(bandwidth), xw(xw), yw(yw), ww(ww), countw(countw), out1(out1), out2(out2),
    kernel(kernel), est(est), flag(flag) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (uword i = begin; i < end; ++i)
    {
      // find the data between out2[i] - bandwidth and out2[i] + bandwidth
      uvec in_windows_1 = linspace<uvec>(0, xw.n_rows-1, xw.n_rows);
      in_windows_1 = in_windows_1.elem(find(all(join_rows(xw.col(1) >= out2(i) - bandwidth(1) - 1e-6,
                                                          xw.col(1) <= out2(i) + bandwidth(1) + 1e-6), 1)));
      if (in_windows_1.n_elem > 0)
      {
        // create model matrix
        mat dx = join_rows(ones<vec>(in_windows_1.n_elem), xw.rows(in_windows_1));
        // compute demean column
        dx.col(2) -= out2(i);
        for (uword j = 0; j < out1.n_elem; ++j)
        {
          // find the data between out1[i] - bandwidth and out1[i] + bandwidth
          uvec range2 = find(all(join_rows(dx.col(1) >= out1(j) - bandwidth(0) - 1e-6,
                                           dx.col(1) <= out1(j) + bandwidth(0) + 1e-6), 1));
          uvec in_windows_2 = in_windows_1.elem(range2);
          if (in_windows_2.n_elem > 0)
          {
            // collect the data in the windwos
            mat lx = dx.rows(range2);
            // compute demean column
            lx.col(1) -= out1(j);
            // collect the data in the windwos
            vec ly = yw(in_windows_2);
            vec lw = ww(in_windows_2);
            vec lc = countw(in_windows_2);
            // find the unique x in the windows and check that the degree of free is sufficient
            mat uni_lx = unique_rows(lx);
            if (uni_lx.n_rows >= 3)
            {
              // compute weights
              vec w = lw % (1-square(lx.col(1)/bandwidth(0))) % (1-square(lx.col(2)/bandwidth(1))) * (9.0/16.0);
              if (kernel == "quar")
                w %= (1 - square(lx.col(1)/bandwidth(0))) % (1 - square(lx.col(2)/bandwidth(1))) * (25.0/16.0);
              /*
                epan: (1 - square(x_minus_range[,1] / bandwidth[1])) *
                      (1 - square(x_minus_range[,2] / bandwidth[2])) * (9.0 / 16.0)
                quar: square(1 - square(x_minus_range[,1] / bandwidth[1])) *
                      square(1 - square(x_minus_range[,2] / bandwidth[2])) * (225.0 / 256.0)
              */
              // fit a WLS
              mat lxw = lx;
              lxw.each_col() %= (w % lc);
              vec p = pinv(lxw.t() * lx) * lx.t() * (w % ly);
              // get the estimation
              est(i, j) = p(0);
            } else
            {
              // return flag = 1 if there is no data in the windows
              flag(i, j) = 1;
            }
          } else
          {
            // return flag = 1 if there is no data in the windows
            flag(i, j) = 1;
          }
        }
      } else
      {
        // return flag = 1 if there is no data in the windows
        flag.row(i).ones();
      }
    }
  }
};

//' Two-dimensional kernel linear local smoother
//'
//' Perform two-dimensional kernel linear local smoother for data \code{(x,y)} with weight \code{w} on \code{xout}.
//'
//' @param bandwidth A numeric vector with two values. The kernel smoothing parameters.
//' @param x A matrix, the variable of of x-axis and y-axis.
//' @param y A vector, the variable of of z-axis. \code{y[i]} is corresponding value of \code{x[i, ]}.
//' @param w A vector, the weight of data. \code{w[i]} is corresponding value of \code{x[i,]}.
//' @param count A vector, the number of observations at \code{x[i, ]}.
//' @param out1 A vector, vector of x-coordinate output time points. It should be a sorted vecotr.
//' @param out2 A vector, vector of y-coordinate output time points. It should be a sorted vecotr.
//' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
//' @return A smoothed covariance estimated by two-dimensional kernel local linear smoother.
//' @examples
//' data("regularExData", package = 'rfda')
//' sparsity <- checkSparsity(regularExData, "sampleID", "t")
//' bwCand <- bwCandChooser(regularExData, "sampleID", "t", sparsity, "gauss", 1)
//' w <- rep(1, nrow(regularExData))
//' bwOpt <- gcvLocPoly1d(bwCand, regularExData$t, regularExData$y, w, "gauss", 0, 1)
//' bwOpt <- rfda:::adjGcvBw1d(bwOpt, sparsity, "gauss", 0)
//' xout <- sort(unique(regularExData$t))
//' meanFunc <- locPoly1d(bwOpt, regularExData$t, regularExData$y, w, xout, "gauss", 0, 1)
//' require(data.table)
//' require(pipeR)
//' demeanDataDT <- merge(data.table(regularExData), data.table(mf = meanFunc, t = xout), by = "t") %>>%
//'   `[`( , `:=`(y = y - mf, variable = "y")) %>>%
//'   setnames(c("t", "y", "sampleID"), c("timePnt", "value", "subId"))
//' RawCov <- rfda:::getRawCrCov(demeanDataDT)
//'
//' xout2 <- seq(min(regularExData$t), max(regularExData$t), len = 30)
//' RawCovNoDiag <- RawCov[t1 != t2]
//' covFunc <- locLinear2d(c(1, 1), as.matrix(RawCovNoDiag[ , .(t1, t2)]), RawCovNoDiag$sse,
//'   RawCovNoDiag$weight, RawCovNoDiag$cnt, xout2, xout2, "gauss")
//' @export
// [[Rcpp::export]]
arma::mat locLinear2d(const arma::vec& bandwidth, const arma::mat& x, const arma::vec& y, const arma::vec& w,
                      const arma::vec& count, const arma::vec& out1, const arma::vec& out2,
                      const std::string& kernel){
  // check data
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  chk_mat(w, "w", "double");
  chk_mat(out1, "out1", "double");
  chk_mat(out2, "out2", "double");
  chk_mat(count, "count", "double");
  if (x.n_rows != y.n_elem || x.n_rows != w.n_elem || x.n_rows != count.n_elem)
    Rcpp::stop("The number of rows of x must be equal to the lengths of y, w and count.\n");
  if (!out1.is_sorted() || !out2.is_sorted())
    Rcpp::stop("out1 and out2 must be sorted.\n");

  // check parameters
  if (kernel != "gauss" && kernel != "gaussvar" && kernel != "epan" && kernel != "quar")
    Rcpp::stop("Unsupported kernel. Kernal must be 'gauss', 'gaussvar', 'epan' or 'quar'.\n");
  if (bandwidth.n_elem != 2 || !is_finite(bandwidth) || any(bandwidth <= 0))
    Rcpp::stop("bandwidth must be a positive numeric vector with 2 elements.\n");

  // remove the data with 0 weights
  uvec actobs = find(w != 0);
  mat xw = x.rows(actobs);
  vec yw = y.elem(actobs), ww = w.elem(actobs), countw = count.elem(actobs);
  // allocate output data
  mat est = zeros<mat>(out2.n_elem, out1.n_elem);
  umat flag = zeros<umat>(out2.n_elem, out1.n_elem);
  // compute parallely the estimates given bandwidth
  if (kernel == "gauss" || kernel == "gaussvar")
  {
    Worker_locLinear2d_gauss locLinear2d(bandwidth, xw, yw, ww, countw, out1, out2, kernel, est);
    RcppParallel::parallelFor(0, out2.n_elem, locLinear2d);
  } else
  {
    Worker_locLinear2d_nongauss locLinear2d(bandwidth, xw, yw, ww, countw, out1, out2, kernel, est, flag);
    RcppParallel::parallelFor(0, out2.n_elem, locLinear2d);
  }
  // return NaN if there is a interval with no data
  mat errOut = {datum::nan};
  if (any(any(flag == 1)))
  {
    RMessage("No enough points in local window, please increase bandwidth.");
    return errOut;
  }
  return est;
}
