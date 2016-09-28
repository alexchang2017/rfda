#include "common.h"
#include "interp.h"
#include <sstream>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

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
      vec x_minus_range = xw - xout(i);
      vec w = ww % exp(-0.5*square(x_minus_range / bandwidth)) / sqrt(2 * datum::pi);
      if (kernel.compare("gaussvar") == 0)
        w %= (1.25 - 0.25 * square(x_minus_range / bandwidth));
      mat dx = ones<mat>(xw.n_elem, degree+1);
      for (uword j = 0; j < degree; ++j)
        dx.col(j+1) = pow(-x_minus_range, j + 1);
      vec p = pinv(dx.t() * (repmat(w, 1, degree+1) % dx)) * dx.t() * (w % yw);
      est(i) = p((uword) drv) * factorial_f(drv) * std::pow(-1.0, drv);
    }
  }
};

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
      uvec in_windows = linspace<uvec>(0, xw.n_elem-1, xw.n_elem);
      in_windows = in_windows.elem(find(all(join_rows(xw <= xout(i) + bandwidth,
                                                      xw >= xout(i) - bandwidth), 1)));
      if (in_windows.n_elem > 0)
      {
        vec lx = xw(in_windows);
        vec ly = yw(in_windows);
        vec lw = ww(in_windows);
        vec uni_lx = unique(lx);
        if (uni_lx.n_elem >= degree + 1)
        {
          vec x_minus_range = lx - xout(i);
          vec w = lw % (1 - square(x_minus_range / bandwidth)) * 0.75;
          if (kernel == "quar")
            w %= (1 - square(x_minus_range / bandwidth)) * (5.0 / 4.0);
          /*
            epan: (1 - square(x_minus_range / bandwidth)) * (3.0 / 4.0)
            quar: square(1 - square(x_minus_range / bandwidth)) * (15.0 / 16.0)
          */
          mat dx = ones<mat>(lx.n_elem, degree+1);
          for (uword j = 0; j < degree; ++j)
            dx.col(j+1) = pow(-x_minus_range, j + 1);
          vec p = pinv(dx.t() * (repmat(w, 1, degree+1) % dx)) * dx.t() * (w % ly);
          est(i) = p((uword) drv) * factorial_f(drv) * std::pow(-1.0, drv);
        } else
        {
          flag(i) = 1;
          est(i) = 0;
        }
      } else
      {
        flag(i) = 1;
        est(i) = 0;
      }
    }
  }
};

//' One-dimensional kernel local polynominal smoother
//'
//' Perform one-dimensional kernel local polynominal smoother for data \code{(x,y)} with weight \code{w} on \code{xout}.
//'
//' @param bandwidth A single numerical value. The kernel smoothing parameter.
//' @param x A vector, the variable of of x-axis.
//' @param y A vector, the variable of of y-axis. \code{y[i]} is corresponding value of \code{x[i]}.
//' @param w A vector, the weight of data. \code{w[i]} is corresponding value of \code{x[i]}.
//' @param xout A vector, vector of output time points. It should be a sorted vecotr.
//' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
//' @param drv An integer, the order of derivative.
//' @param degree An integer, the degree of polynomial.
//' @return A estimated value on \code{xout} by one-dimensional kernel local polynominal smoother.
//' @examples
//' x <- runif(100, 0, 10)
//' y <- rnorm(100)
//' xout <- sort(runif(200, 0, 10))
//' est <- locPoly1d(1.2, x, y, rep(1, 100), xout, 'gauss', 0, 1)
//' require(ggplot2)
//' ggplot(data.frame(x,y), aes(x,y)) + geom_point() +
//'   geom_line(aes(xout, est), data = data.frame(xout, est))
//' @export
// [[Rcpp::export]]
arma::vec locPoly1d(const double& bandwidth, const arma::vec& x, const arma::vec& y,
                    const arma::vec& w, const arma::vec& xout, const std::string& kernel,
                    const double& drv, const double& degree){
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  chk_mat(w, "w", "double");
  if (x.n_elem != y.n_elem || x.n_elem != w.n_elem)
    Rcpp::stop("The lengths of x, y and w must be equal.\n");
  chk_mat(xout, "xout", "double");
  if (!xout.is_sorted())
    Rcpp::stop("xout must be sorted.\n");

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

  uvec actobs = find(w != 0);
  vec xw = x.elem(actobs), yw = y.elem(actobs), ww = w.elem(actobs), est = zeros<vec>(xout.n_elem);
  uvec flag = zeros<uvec>(xout.n_elem);
  if (kernel == "gauss" || kernel == "gaussvar")
  {
    Worker_locPoly1d_gauss locPoly1d(bandwidth, xw, yw, ww, xout, kernel, drv, degree, est);
    RcppParallel::parallelFor(0, xout.n_elem, locPoly1d);
  } else
  {
    Worker_locPoly1d_nongauss locPoly1d(bandwidth, xw, yw, ww, xout, kernel, drv, degree, est, flag);
    RcppParallel::parallelFor(0, xout.n_elem, locPoly1d);
  }
  vec errOut = {datum::nan};
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
//' @param bwCand A numerical vector for the candidates of bandwidth
//' @param x A vector, the variable of of x-axis.
//' @param y A vector, the variable of of y-axis. \code{y[i]} is corresponding value of \code{x[i]}.
//' @param w A vector, the weight of data. \code{w[i]} is corresponding value of \code{x[i]}.
//' @param kernel A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.
//' @param drv An integer, the order of derivative.
//' @param degree An integer, the degree of polynomial.
//' @return A optimal bandwidth selected by minimizing gcv scores.
//' @examples
//' data("regularExData", package = 'rfda')
//' regBwCand <- bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1)
//' w <- rep(1, nrow(regularExData))
//' bw_opt <- gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y, w, "gauss", 0, 1)
//' @export
// [[Rcpp::export]]
double gcv_locPoly1d(arma::vec bwCand, const arma::vec& x, const arma::vec& y,
                     const arma::vec& w, const std::string& kernel,
                     const double& drv, const double& degree){
  chk_mat(bwCand, "bwCand", "double");
  chk_mat(x, "x", "double");
  chk_mat(y, "y", "double");
  chk_mat(w, "w", "double");
  if (x.n_elem != y.n_elem || x.n_elem != w.n_elem)
    Rcpp::stop("The lengths of x, y and w must be equal.\n");

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

  vec k0 = kernelDensity(zeros<vec>(1.0), kernel);
  double r = x.max() - x.min(), bw_opt;
  vec gcv_param = pow(1 - r * k0(0) / bwCand / (double) x.n_elem, 2.0),
    gcv = zeros<vec>(bwCand.n_elem), new_est(x.n_elem);

  bool con = true, sparse = false, secondRun = false;
  std::string interp1_method = "spline";
  vec xout = sort(unique(x)), xout2 = xout;
  if (xout.n_elem > 101)
    xout2 = linspace<vec>(min(xout), max(xout), 101);
  vec est = zeros<vec>(xout2.n_elem);
  while (con)
  {
    gcv.fill(1e35);
    for (uword k = 0; k < bwCand.n_elem; ++k)
    {
      est = locPoly1d(bwCand(k), x, y, w, xout2, kernel, drv, degree);
      if (is_finite(est)){
        if (xout2.n_elem == xout.n_elem){
          new_est.zeros();
          for (uword k = 0; k < xout.n_elem; ++k)
            new_est.elem(find(x == xout(k))).fill(est(k));
        } else {
          new_est = interp1(xout2, est, x, interp1_method);
        }
      }

      if (is_finite(new_est))
      {
        gcv(k) = dot(y - new_est, y - new_est) / gcv_param(k);
        if (k > 0 && gcv(k) > gcv(k-1))
        {
          con = false;
          break;
        }
      }
    }

    if (all(gcv == 1e35))
    {
      if (!secondRun && bwCand(9) < r)
      {
        bw_opt = bwCand(9);
        sparse = true;
      } else
      {
        Rcpp::stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n");
      }
    } else
    {
      uword min_gcv_idx = as_scalar(find(gcv == min(gcv), 1, "first"));
      bw_opt = bwCand(min_gcv_idx);
    }

    if (bw_opt == r)
    {
      con = false;
      Rcpp::stop("The data is too sparse, optimal bandwidth includes all the data! You may want to change to Gaussian kernel!\n");
    } else if (bw_opt == bwCand(9) && !secondRun)
    {
      uvec gcvInf = find(gcv == 1e35);
      double minBW;
      if (sparse || gcvInf.n_elem == 9)
      {
        RMessage("The data is too sparse, retry with larger bandwidths!");
        minBW = bwCand(9) * 1.01;
      } else
      {
        RMessage("Bandwidth candidates are too small, retry with larger choices now!");
        minBW = bwCand(8);
      }
      vec newr = r * linspace<vec>(0.5, 1, 11);
      uword min_idx = as_scalar(find(minBW < newr, 1, "first"));
      double q = std::pow(newr(min_idx) / minBW, 1.0 / 9.0);
      bwCand = linspace<vec>(0, 9, 10);
      for (uword i = 0; i < bwCand.n_elem; ++i)
        bwCand(i) = std::pow(q, bwCand(i)) * minBW;
      bwCand = sort(bwCand);

      RMessage("New bwmu candidates: ");
      std::stringstream bwCand_str;
      bwCand.print(bwCand_str);
      RMessage(bwCand_str.str());
    } else if (bw_opt < bwCand(9) || secondRun)
    {
      con = false;
    }
    secondRun = true;
  }

  return bw_opt;
}

