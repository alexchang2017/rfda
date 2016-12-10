#include "common.h"
#include "interp.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct Worker_locLinearRotate2d_gauss: public RcppParallel::Worker {
  const vec& bandwidth;
  const mat& xw;
  const vec& yw;
  const vec& ww;
  const vec& countw;
  const  mat& outMat;
  const std::string& kernel;
  vec& est;

  Worker_locLinearRotate2d_gauss(const vec& bandwidth, const mat& xw, const vec& yw, const vec& ww, const vec& countw,
                                 const mat& outMat, const std::string& kernel, vec& est):
    bandwidth(bandwidth), xw(xw), yw(yw), ww(ww), countw(countw), outMat(outMat), kernel(kernel), est(est) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (uword i = begin; i < end; ++i)
    {
      // create model matrix
      mat dx = join_rows(ones<vec>(xw.n_rows), xw);
      // compute demean columns
      dx.col(1) -= outMat(i, 0);
      dx.col(2) -= outMat(i, 1);
      vec w = ww % exp(-0.5 * square(dx.col(1) / bandwidth(0))) %
        exp(-0.5 * square(dx.col(2) / bandwidth(1))) / (2 * datum::pi);
      if (kernel == "gaussvar")
        w %= (1.25 - 0.25 * square(dx.col(1) / bandwidth(0))) %
          (1.5 - 0.5 * square(dx.col(2) / bandwidth(1)));
      dx.col(1) = square(dx.col(1));
      mat dxw = dx;
      dxw.each_col() %= w;
      vec p = pinv(dxw.t() * dx) * dx.t() * (w % yw / countw);
      est(i) = p(0);
    }
  }
};

struct Worker_locLinearRotate2d_nongauss: public RcppParallel::Worker {
  const vec& bandwidth;
  const mat& xw;
  const vec& yw;
  const vec& ww;
  const vec& countw;
  const mat& outMat;
  const std::string& kernel;
  vec& est;
  uvec& flag;

  Worker_locLinearRotate2d_nongauss(const vec& bandwidth, const mat& xw, const vec& yw, const vec& ww,
                              const vec& countw, const mat& outMat, const std::string& kernel,
                              vec& est, uvec& flag):
    bandwidth(bandwidth), xw(xw), yw(yw), ww(ww), countw(countw), outMat(outMat),
    kernel(kernel), est(est), flag(flag) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (uword i = begin; i < end; ++i)
    {
      // find the data between out2[i] - bandwidth and out2[i] + bandwidth
      uvec in_windows = linspace<uvec>(0, xw.n_rows-1, xw.n_rows);
      in_windows = in_windows.elem(find(all(join_rows(
        join_rows(xw.col(0) >= outMat(i, 0) - bandwidth(1) - 1e-6, xw.col(0) <= outMat(i, 0) + bandwidth(1) + 1e-6),
        join_rows(xw.col(1) >= outMat(i, 1) - bandwidth(1) - 1e-6, xw.col(1) <= outMat(i, 1) + bandwidth(1) + 1e-6)), 1)));
      if (in_windows.n_elem >= 2)
      {
        // create model matrix
        mat lx = join_rows(ones<vec>(in_windows.n_elem), xw.rows(in_windows));
        // compute demean columns
        lx.col(1) -= outMat(i, 0);
        lx.col(2) -= outMat(i, 1);
        vec ly = yw(in_windows);
        vec lw = ww(in_windows);
        vec lc = countw(in_windows);

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
        lx.col(1) = square(lx.col(1));
        mat lxw = lx;
        lxw.each_col() %= w;
        vec p = pinv(lxw.t() * lx) * lx.t() * (w % ly / lc);
        // get the estimation
        est(i) = p(0);
      } else if (in_windows.n_elem == 1)
      {
        est(i) = as_scalar(yw(in_windows));
      } else
      {
        // return flag = 1 if there is no data in the windows
        flag(i) = 1;
      }
    }
  }
};

// [[Rcpp::export]]
arma::vec locLinearRotate2d_cpp(const arma::vec& bandwidth, const arma::mat& x, const arma::vec& y, const arma::vec& w,
                                const arma::vec& count, const arma::mat outMat, const std::string& kernel){
  // check data
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  chk_mat(w, "w", "double");
  chk_mat(outMat, "outMat", "double");
  chk_mat(count, "count", "double");
  if (x.n_rows != y.n_elem || x.n_rows != w.n_elem || x.n_rows != count.n_elem)
    Rcpp::stop("The number of rows of x must be equal to the lengths of y, w and count.\n");
  if (x.n_cols != 2)
    Rcpp::stop("x must be a matrix with 2 columns.\n");
  if (outMat.n_cols != 2)
    Rcpp::stop("outMat must be a matrix with 2 columns.\n");

  // check parameters
  if (kernel != "gauss" && kernel != "gaussvar" && kernel != "epan" && kernel != "quar")
    Rcpp::stop("Unsupported kernel. Kernal must be 'gauss', 'gaussvar', 'epan' or 'quar'.\n");
  if (bandwidth.n_elem != 2 || !is_finite(bandwidth) || any(bandwidth <= 0))
    Rcpp::stop("bandwidth must be a positive numeric vector with 2 elements.\n");

  // rotate matrix
  mat R = std::sqrt(2)/2 * ones<mat>(2, 2);
  R(1, 0) *= -1.0;

  // remove the data with 0 weights
  uvec actobs = find(w != 0);
  mat xw = x.rows(actobs) * R, outMat2 = outMat * R;
  vec yw = y.elem(actobs), ww = w.elem(actobs), countw = count.elem(actobs);
  // allocate output data
  vec est = zeros<vec>(outMat2.n_rows);
  uvec flag = zeros<uvec>(outMat2.n_rows);
  // compute parallely the estimates given bandwidth
  if (kernel == "gauss" || kernel == "gaussvar")
  {
    Worker_locLinearRotate2d_gauss locLinearRotate2d(bandwidth, xw, yw, ww, countw, outMat2, kernel, est);
    RcppParallel::parallelFor(0, outMat2.n_rows, locLinearRotate2d);
  } else
  {
    Worker_locLinearRotate2d_nongauss locLinearRotate2d(bandwidth, xw, yw, ww, countw, outMat2, kernel, est, flag);
    RcppParallel::parallelFor(0, outMat2.n_rows, locLinearRotate2d);
  }
  // return NaN if there is a interval with no data
  vec errOut(1);
  errOut(0) = datum::nan;
  if (any(flag == 1))
  {
    RMessage("No enough points in local window, please increase bandwidth.");
    return errOut;
  }
  return est;
}
