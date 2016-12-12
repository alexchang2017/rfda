#include "common.h"
#include "interp.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

// a parallel worker to compute estimate of local kernel polynominal smoother with gaussian-like kernel on (x, y)
struct Worker_locPoly1d_gauss: public RcppParallel::Worker {
  const double& bandwidth;
  const vec& xw;
  const vec& yw;
  const vec& ww;
  const vec& xout;
  const std::string& kernel;
  const double& drv;
  const double& degree;
  vec& est;

  Worker_locPoly1d_gauss(const double& bandwidth, const vec& xw, const vec& yw, const vec& ww,
                         const vec& xout, const std::string& kernel, const double& drv,
                         const double& degree, vec& est):
    bandwidth(bandwidth), xw(xw), yw(yw), ww(ww), xout(xout), kernel(kernel),
    drv(drv), degree(degree), est(est) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (uword i = begin; i < end; ++i)
    {
      // compute demean data
      vec x_minus_range = xw - xout(i);
      // compute weights
      vec w = ww % exp(-0.5*square(x_minus_range / bandwidth)) / sqrt(2 * datum::pi);
      if (kernel == "gaussvar")
        w %= (1.25 - 0.25 * square(x_minus_range / bandwidth));
      // get model matrix
      mat dx = ones<mat>(xw.n_elem, degree+1);
      for (uword j = 0; j < degree; ++j)
        dx.col(j+1) = pow(-x_minus_range, j + 1);
      // fit a WLS
      vec p = pinv(dx.t() * (repmat(w, 1, degree+1) % dx)) * dx.t() * (w % yw);
      // get the estimation
      est(i) = p((uword) drv) * factorial_f(drv) * std::pow(-1.0, drv);
    }
  }
};

// a parallel worker to compute estimate of local kernel polynominal smoother with non-gaussian kernel on (x, y)
struct Worker_locPoly1d_nongauss: public RcppParallel::Worker {
  const double& bandwidth;
  const vec& xw;
  const vec& yw;
  const vec& ww;
  const vec& xout;
  const std::string& kernel;
  const double& drv;
  const double& degree;
  vec& est;
  uvec& flag;

  Worker_locPoly1d_nongauss(const double& bandwidth, const vec& xw, const vec& yw, const vec& ww,
                         const vec& xout, const std::string& kernel, const double& drv,
                         const double& degree, vec& est, uvec& flag):
    bandwidth(bandwidth), xw(xw), yw(yw), ww(ww), xout(xout), kernel(kernel),
    drv(drv), degree(degree), est(est), flag(flag) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (uword i = begin; i < end; ++i)
    {
      // find the data between xout[i] - bandwidth and xout[i] + bandwidth
      uvec in_windows = linspace<uvec>(0, xw.n_elem-1, xw.n_elem);
      in_windows = in_windows.elem(find(all(join_rows(xw <= xout(i) + bandwidth,
                                                      xw >= xout(i) - bandwidth), 1)));
      if (in_windows.n_elem > 0)
      {
        // collect the data in the windwos
        vec lx = xw.elem(in_windows);
        vec ly = yw.elem(in_windows);
        vec lw = ww.elem(in_windows);
        // find the unique x in the windows and check that the degree of free is sufficient
        vec uni_lx = unique(lx);
        if (uni_lx.n_elem >= degree + 1)
        {
          // compute demean data
          vec x_minus_range = lx - xout(i);
          // compute weights
          vec w = lw % (1 - square(x_minus_range / bandwidth)) * 0.75;
          if (kernel == "quar")
            w %= (1 - square(x_minus_range / bandwidth)) * (5.0 / 4.0);
          /*
            epan: (1 - square(x_minus_range / bandwidth)) * (3.0 / 4.0)
            quar: square(1 - square(x_minus_range / bandwidth)) * (15.0 / 16.0)
          */
          // get model matrix
          mat dx = ones<mat>(lx.n_elem, degree+1);
          for (uword j = 0; j < degree; ++j)
            dx.col(j+1) = pow(-x_minus_range, j + 1);
          // fit a WLS
          vec p = pinv(dx.t() * (repmat(w, 1, degree+1) % dx)) * dx.t() * (w % ly);
          // get the estimation
          est(i) = p((uword) drv) * factorial_f(drv) * std::pow(-1.0, drv);
        } else
        {
          // return flag = 1 if there is no data in the windows
          flag(i) = 1;
        }
      } else
      {
        // return flag = 1 if there is no data in the windows
        flag(i) = 1;
      }
    }
  }
};

// [[Rcpp::export]]
arma::vec locPoly1d_cpp(const double& bandwidth, const arma::vec& x, const arma::vec& y,
                        const arma::vec& w, const arma::vec& xout, const std::string& kernel,
                        const double& drv, const double& degree){
  // check data
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  chk_mat(w, "w", "double");
  chk_mat(xout, "xout", "double");
  if (x.n_elem != y.n_elem || x.n_elem != w.n_elem)
    Rcpp::stop("The lengths of x, y and w must be equal.\n");
  if (!xout.is_sorted())
    Rcpp::stop("xout must be sorted.\n");

  // check parameters
  if (kernel != "gauss" && kernel != "gaussvar" && kernel != "epan" && kernel != "quar")
    Rcpp::stop("Unsupported kernel. Kernal must be 'gauss', 'gaussvar', 'epan' or 'quar'.\n");
  if (!is_finite(bandwidth) || bandwidth <= 0)
    Rcpp::stop("bandwidth must be a positive number.\n");
  if (!is_finite(degree))
    Rcpp::stop("degree must be a nonnegative number.\n");
  if (degree < 0 || std::abs(degree - floor(degree)) > 1e-6)
    Rcpp::stop("degree must be a nonnegative number.\n");
  if (!is_finite(drv))
    Rcpp::stop("drv must be a nonnegative number.\n");
  if (drv < 0 || std::abs(drv - floor(drv)) > 1e-6)
    Rcpp::stop("drv must be a nonnegative number.\n");
  if (degree < drv)
    Rcpp::stop("Degree of Polynomial should be not less than the order of derivative.\n");

  // remove the data with 0 weights
  uvec actObs = find(w != 0);
  vec xw = x.elem(actObs), yw = y.elem(actObs), ww = w.elem(actObs), est = zeros<vec>(xout.n_elem);
  // allocate output data
  uvec flag = zeros<uvec>(xout.n_elem);
  // compute parallely the estimates given bandwidth
  if (kernel == "gauss" || kernel == "gaussvar")
  {
    Worker_locPoly1d_gauss locPoly1d_worker(bandwidth, xw, yw, ww, xout, kernel, drv, degree, est);
    RcppParallel::parallelFor(0, xout.n_elem, locPoly1d_worker);
  } else
  {
    Worker_locPoly1d_nongauss locPoly1d_worker(bandwidth, xw, yw, ww, xout, kernel, drv, degree, est, flag);
    RcppParallel::parallelFor(0, xout.n_elem, locPoly1d_worker);
  }
  // return NaN if there is a interval with no data
  vec errOut(1);
  errOut(0) = datum::nan;
  if (any(flag == 1))
  {
    RMessage("Too many gaps, please increase bandwidth.");
    return errOut;
  }
  return est;
}

//' Find the optimal bandwidth for one-dimensional kernel local polynominal smoother
//'
//' Find the optimal bandwidth used in \code{\link{locPoly1d}}.
//'
//' @param bwCand A numerical vector for the candidates of bandwidth.
//' @param x A vector, the variable of of x-axis.
//' @param y A vector, the variable of of y-axis. \code{y[i]} is corresponding value of \code{x[i]}.
//' @param w A vector, the weight of data. \code{w[i]} is corresponding value of \code{x[i]}.
//' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
//' @param drv An integer, the order of derivative.
//' @param degree An integer, the degree of polynomial.
//' @return A optimal bandwidth selected by minimizing gcv scores.
//' @examples
//' data("regularExData", package = 'rfda')
//' bwCand <- bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1)
//' w <- rep(1, nrow(regularExData))
//' bwOpt <- gcvLocPoly1d(bwCand, regularExData$t, regularExData$y, w, "gauss", 0, 1)
//' bwOpt <- adjGcvBw(bwOpt, 2, "gauss", 0)
//' @export
// [[Rcpp::export]]
double gcvLocPoly1d(arma::vec bwCand, const arma::vec& x, const arma::vec& y,
                    const arma::vec& w, const std::string& kernel,
                    const double& drv, const double& degree){
  // check data
  chk_mat(bwCand, "bwCand", "double");
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  chk_mat(w, "w", "double");
  if (any(bwCand <= 0))
    Rcpp::stop("The elements in bwCand must be greater 0.\n");
  if (x.n_elem != y.n_elem || x.n_elem != w.n_elem)
    Rcpp::stop("The lengths of x, y and w must be equal.\n");

  // check parameters
  if (kernel != "gauss" && kernel != "gaussvar" && kernel != "epan" && kernel != "quar")
    Rcpp::stop("Unsupported kernel. Kernal must be 'gauss', 'gaussvar', 'epan' or 'quar'.\n");
  if (!is_finite(degree))
    Rcpp::stop("degree must be a nonnegative number.\n");
  if (degree < 0 || std::abs(degree - floor(degree)) > 1e-6)
    Rcpp::stop("degree must be a nonnegative number.\n");
  if (!is_finite(drv))
    Rcpp::stop("drv must be a nonnegative number.\n");
  if (drv < 0 || std::abs(drv - floor(drv)) > 1e-6)
    Rcpp::stop("drv must be a nonnegative number.\n");
  if (degree < drv)
    Rcpp::stop("Degree of Polynomial should be not less than the order of derivative.\n");

  // find the GCV-adjusted term
  vec k0 = kernelDensity(zeros<vec>(1), kernel);
  // find the range of x and allocate output
  double r = x.max() - x.min(), bwOpt;
  // find the GCV-adjusted term
  vec gcv_param = pow(1 - r * k0(0) / bwCand / (double) x.n_elem, 2.0);

  // initialize gcv
  vec gcv = zeros<vec>(bwCand.n_elem);
  bool con = true, sparse = false, secondRun = false;
  std::string interp1_method = "spline";
  vec xout = sort(unique(x)), xout2 = xout;
  if (xout.n_elem > 101)
    xout2 = linspace<vec>(min(xout), max(xout), 101);
  vec est = zeros<vec>(xout2.n_elem), new_est = zeros<vec>(x.n_elem), diff = zeros<vec>(x.n_elem);
  uvec nonfiniteLoc;
  while (con)
  {
    // fill gcv score with big number
    gcv.fill(datum::inf);
    for (uword k = 0; k < bwCand.n_elem; ++k)
    {
      // get smoothed values
      est = locPoly1d_cpp(bwCand(k), x, y, w, xout2, kernel, drv, degree);
      if (is_finite(est)){
        // map the smoothed values onto x
        if (xout2.n_elem == xout.n_elem){
          new_est.zeros();
          for (uword k = 0; k < xout.n_elem; ++k)
            new_est.elem(find(x == xout(k))).fill(est(k));
        } else {
          new_est = interp1_cpp(xout2, est, x, interp1_method);
        }
      } else {
        new_est.fill(datum::nan);
      }

      if (is_finite(new_est))
      {
        // compute gcv scores
        diff = y - new_est;
        gcv(k) = dot(diff, diff) / gcv_param(k);
        if (k > 0 && gcv(k) > gcv(k-1))
        {
          con = false;
          break;
        }
      }
    }

    nonfiniteLoc = find_nonfinite(gcv);
    if (nonfiniteLoc.n_elem == bwCand.n_elem)
    {
      // if bw are all too small, then re-compute or
      // stop implementation if it has the same situation in second round
      if (!secondRun && bwCand(bwCand.n_elem - 1) < r)
      {
        bwOpt = bwCand(bwCand.n_elem - 1);
        sparse = true;
      } else
      {
        Rcpp::stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n");
      }
    } else
    {
      // select the optimal bandwidth with minimum gcv scores
      uword min_gcv_idx = as_scalar(find(gcv == min(gcv), 1, "first"));
      bwOpt = bwCand(min_gcv_idx);
    }

    if (bwOpt == r)
    {
      // stop implementation if data is too sparse
      con = false;
      Rcpp::stop("The data is too sparse, optimal bandwidth includes all the data! You may want to change to Gaussian kernel!\n");
    } else if (bwOpt == bwCand(bwCand.n_elem - 1) && !secondRun)
    {
      // re-compute with new bandwidth candidates
      double minBW;
      if (sparse || nonfiniteLoc.n_elem == bwCand.n_elem-1)
      {
        RMessage("The data is too sparse, retry with larger bandwidths!");
        minBW = bwCand(bwCand.n_elem - 1) * 1.01;
      } else
      {
        RMessage("Bandwidth candidates are too small, retry with larger choices now!");
        minBW = bwCand(bwCand.n_elem - 2);
      }
      vec newr = r * linspace<vec>(0.5, 1, bwCand.n_elem + 1);
      uword min_idx = as_scalar(find(minBW < newr, 1, "first"));
      double q = std::pow(newr(min_idx) / minBW, 1.0 / ((double) bwCand.n_elem - 1.0));
      bwCand = linspace<vec>(0, bwCand.n_elem - 1, bwCand.n_elem);
      for (uword i = 0; i < bwCand.n_elem; ++i)
        bwCand(i) = std::pow(q, bwCand(i)) * minBW;
      bwCand = sort(bwCand);

      RMessage("New bwmu candidates: ");
      std::stringstream bwCand_str;
      bwCand.print(bwCand_str);
      RMessage(bwCand_str.str());
    } else if (bwOpt < bwCand(bwCand.n_elem - 1) || secondRun)
    {
      con = false;
    }
    secondRun = true;
  }
  return bwOpt;
}

