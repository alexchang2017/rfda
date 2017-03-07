#include "common.h"
#include <algorithm>

arma::umat lookup(const arma::vec& edges, const arma::mat& x) {
  if (!edges.is_sorted())
    Rcpp::stop("edges is not strictly monotonic increasing.");

  umat idx(size(x));
  const double* pos;
  for (uword i = 0; i < x.n_elem; ++i) {
    pos = std::upper_bound(edges.begin(), edges.end(), x(i));
    idx(i) = std::distance(edges.begin(), pos);
  }
  return idx;
}

// spline function in C++ implementation
// [[Rcpp::export]]
arma::mat spline_cpp(const arma::vec& x, const arma::mat& y, const arma::vec& xi){
  chk_mat(x, "x");
  chk_mat(y, "y");
  chk_mat(xi, "xi");

  uword n = x.n_elem;
  if (n < 2)
    Rcpp::stop("spline: requires at least 2 points.");

  mat a = y;
  uvec szy;
  szy << y.n_rows << y.n_cols << endr;
  if (szy(1) == n && szy(0) != n+2 && szy(0) != n && szy(0) >= 1)
    inplace_trans(a);
  if (szy(1) == n+2 && szy(0) != n+2 && szy(0) != n && szy(0) >= 1)
    inplace_trans(a);
  if (x.n_elem != a.n_rows && a.n_rows != x.n_elem+2)
    Rcpp::stop("The number of rows of y must be equal to the length of x.\n");

  bool complete = false;
  rowvec dfs, dfe;
  if (a.n_rows == n+2) {
    complete = true;
    dfs = a.row(0);
    dfe = a.row(a.n_rows-1);
    a = a.rows(1, a.n_rows-2);
  }

  if (!x.is_sorted())
    RMessage("x are not strictly monotonic increasing.\nx will be sorted.");
  uvec xu_idx = find_unique(x);
  vec xu = x(xu_idx);
  mat au = a.rows(xu_idx);
  if (xu.n_elem != x.n_elem) {
    RMessage("The grid vectors are not strictly monotonic increasing.");
    RMessage("The values of y for duplicated values of x will be averaged.");
    for (uword k = 0; k < xu.n_elem; ++k)
      au.row(k) = mean(a.rows(find(x == xu(k))));
    n = xu_idx.n_elem;
  }
  if (!xu.is_sorted()) {
    uvec si = sort_index(xu);
    xu = xu(si);
    au = au.rows(si);
  }

  mat ca, cb, cc, cd;
  vec xou = xu, h = diff(xu);
  if (complete) {
    if (n == 2) {
      cd = (dfs + dfe) / std::pow(xu(1) - xu(0), 2.0) +
        2.0 * (au.row(0) - au.row(1)) / std::pow(xu(1) - xu(0), 3.0);
      cc = (-2.0 * dfs - dfe) / (xu(1) - xu(0)) -
        3.0 * (au.row(0) - au.row(1)) / std::pow(xu(1) - xu(0), 2.0);
      cb = dfs;
      ca = au.row(0);
    } else {
      mat g = zeros<mat>(n, au.n_cols);
      g.row(0) = (au.row(1) - au.row(0)) / h(0) - dfs;
      g.rows(1, n-2) = (au.rows(2, n-1) - au.rows(1, n-2)) / repmat(h.subvec(1, n-2), 1, au.n_cols) -
        (au.rows(1, n-2) - a.rows(0, n-3)) / repmat(h.subvec(0, n-3), 1, au.n_cols);
      g.row(n-1) = dfe-(au.row(n-1) - au.row(n-2)) / h(n-2);

      ca = au;
      cc = solve(diagmat(h/6.0, -1) +
        diagmat(join_cols(join_cols(h(0)/3 * ones<vec>(1), (h.head(n-2) + h.tail(n-2))/3),
                          h(n-2)/3.0*ones<vec>(1))) + diagmat(h/6.0, 1), 0.5 * g);
      cb = diff(au) / repmat(h.head(n-1), 1, au.n_cols) -
        repmat(h.head(n-1), 1, au.n_cols) / 3.0 % (cc.rows(1, n-1) + 2 * cc.rows(0, n-2));
      cd = diff(cc) / (3.0 * repmat(h.head(n-1), 1, au.n_cols));
      ca = ca.head_rows(n-1);
      cb = cb.head_rows(n-1);
      cc = cc.head_rows(n-1);
      cd = cd.head_rows(n-1);
    }
  } else {
    if (n == 2) {
      cd.zeros(1, au.n_cols);
      cc.zeros(1, au.n_cols);
      cb = (au.row(1) - au.row(0)) / (xu(1) - xu(0));
      ca = au.row(0);
    } else if (n == 3) {
      n = 2;
      cd.zeros(1, au.n_cols);
      cc = (au.row(0) - au.row(2)) / ((xu(2) - xu(0)) * (xu(1) - xu(2))) +
        (au.row(1) - au.row(0)) / ((xu(1) - xu(0)) * (xu(1) - xu(2)));
      cb = (au.row(1) - au.row(0)) * (xu(2) - xu(0)) /  ((xu(1) - xu(0)) * (xu(2) - xu(1))) +
        (au.row(0) - au.row(2)) * (xu(1) - xu(0)) /  ((xu(2) - xu(0)) * (xu(2) - xu(1)));
      ca = au.row(0);
      xou << min(x) << max(x) << endr;
    } else {
      mat g = zeros<mat>(n-2, au.n_cols);
      g.row(0) = 3.0 / (h(0) + h(1)) *
        (au.row(2) - au.row(1) - h(1) / h(0) * (au.row(1) - au.row(0)));
      g.row(n-3) = 3.0 / (h(n-2) + h(n-3)) *
        (h(n-3) / h(n-2) * (au.row(n-1) - au.row(n-2)) - (au.row(n-2) - au.row(n-3)));

      if (n > 4) {
        cc.zeros(n, au.n_cols);
        g.rows(1, n-4) = 3.0 * diff(au.rows(2, n-2)) / repmat(h.subvec(2, n-3), 1, au.n_cols) -
          3.0 * diff(au.rows(1, n-3)) / repmat(h.subvec(1, n-4), 1, au.n_cols);

        vec dg = 2.0 * (h.head(n-2) + h.tail(n-2)),
          ldg = h.subvec(1, n-3), udg = h.subvec(1, n-3);
        dg(0) = dg(0) - h(0);
        dg(n-3) = dg(n-3) - h(n-2);
        udg(0) = udg(0) - h(0);
        ldg(n-4) = ldg(n-4) - h(n-2);
        cc.rows(1, n-2) = solve(diagmat(ldg, -1) + diagmat(dg) + diagmat(udg, 1), g);
      } else {
        cc.zeros(n, au.n_cols);
        mat tmp;
        tmp << h(0) + 2.0 * h(1) << h(1) - h(0) << endr
            << h(1) - h(2) << 2.0 * h(1) + h(2) << endr;
        cc.rows(1, 2) = solve(tmp, g);
      }

      ca = au;
      cc.row(0) = cc.row(1) + h(0) / h(1) * (cc.row(1) - cc.row(2));
      cc.row(n-1) = cc.row(n-2) + h(n-2) / h(n-3) * (cc.row(n-2) - cc.row(n-3));
      cb = diff(ca);
      cb.each_col() /= h.head(n-1);
      cb -= repmat(h.head(n-1), 1, au.n_cols) / 3.0 % (cc.rows(1, n-1) + 2.0 * cc.rows(0, n-2));
      cd = diff(cc) / 3.0;
      cd.each_col() /= h.head(n-1);
      ca = ca.head_rows(n-1);
      cb = cb.head_rows(n-1);
      cc = cc.head_rows(n-1);
      cd = cd.head_rows(n-1);
    }
  }

  uvec idx = zeros<uvec>(xi.n_elem);
  for (uword i = 1; i < xou.n_elem-1; ++i)
    idx.elem(find(xou(i) <= xi)).fill(i);
  mat s_mat = repmat(xi - xou.elem(idx), 1, au.n_cols);
  mat ret = ca.rows(idx) + s_mat % cb.rows(idx) + square(s_mat) % cc.rows(idx) +
    pow(s_mat, 3) % cd.rows(idx);
  return ret;
}

// interp1 function in C++ implementation
// [[Rcpp::export]]
arma::mat interp1_cpp(const arma::vec& x, const arma::mat& y, const arma::vec& xi,
                      const std::string& method){
  chk_mat(x, "x");
  chk_mat(y, "y");
  chk_mat(xi, "xi");

  uword n = x.n_elem;
  mat a = y;
  uvec szy;
  szy << y.n_rows << y.n_cols << endr;
  if (szy(1) == n && szy(0) != n+2 && szy(0) != n && szy(0) != 1)
    inplace_trans(a);
  if (x.n_elem != a.n_rows)
    Rcpp::stop("The number of rows of y must be equal to the length of x.\n");

  if (!x.is_sorted())
    RMessage("x are not strictly monotonic increasing.\nx will be sorted.");
  uvec xu_idx = find_unique(x);
  vec xu = x(xu_idx);
  mat au = a.rows(xu_idx);
  if (xu.n_elem != x.n_elem) {
    RMessage("The grid vectors are not strictly monotonic increasing.");
    RMessage("The values of y for duplicated values of x will be averaged.");
    for (uword k = 0; k < xu.n_elem; ++k)
      au.row(k) = mean(a.rows(find(x == xu(k))));
  }
  if (!xu.is_sorted()) {
    uvec si = sort_index(xu);
    xu = xu(si);
    au = au.rows(si);
  }

  mat yi;
  if (method.compare("linear") == 0) {
    if (x.n_elem <= 1)
      Rcpp::stop("interp1 - linear: requires at least two non-NA values.\n");
    mat cb = diff(au) / repmat(diff(xu), 1, au.n_cols);
    mat ca = au.rows(0, xu.n_elem-2);
    uvec idx = zeros<uvec>(xi.n_elem);
    for (uword i = 1; i < xu.n_elem-1; ++i)
      idx.elem(find(xu(i) <= xi)).fill(i);
    mat s_mat = repmat(xi - xu.elem(idx), 1, au.n_cols);
    yi = ca.rows(idx) + s_mat % cb.rows(idx);
  } else if (method.compare("spline") == 0) {
    yi = spline_cpp(xu, au, xi);
  } else {
    Rcpp::stop("Method only support linear and spline.\n");
  }
  return yi;
}

//' 2-D data interpolation.
//'
//' Returns interpolated values of a 2-D function at specific query points using
//' linear interpolation. The extrapolation is used, please be caution in using the
//' values which \code{xi} is larger than \code{max(x)/max(y)} and smaller than \code{min(x)/min(y)}.
//'
//' @param x A vector with n1 elements, \code{x[i]} is a support, \code{i = 1, ..., n1}.
//'   If \code{x} is not sorted, it will be sorted. If \code{x} is not unique, the corresponding \code{v} values
//'   will be averaged.
//' @param y A vector with n2 elements, \code{y[j]} is a support, \code{j = 1, ..., n2}.
//'   If \code{y} is not sorted, it will be sorted. If \code{y} is not unique, the corresponding \code{v} values
//'   will be averaged.
//' @param v A matrix with size n1 by n2, \code{v[i, j]} is the corresponding value at grid \code{(x[i], y[j])}.
//' @param xi A vector with m elements, \code{xi[k]} is the point which you want to interpolate,
//'   \code{k = 1, ..., m1}.
//' @param yi A vector with m elements, \code{yi[l]} is the point which you want to interpolate,
//'   \code{l = 1, ..., m2}.
//' @param method A string "linear" or "spline", the method of interpolation.
//' @return A matrix with the interpolated values corresponding to \code{xi} and \code{yi}.
//' @section Reference:
//' Cleve Moler, Numerical Computing with MATLAB, chapter 3,
//'   \url{http://www.mathworks.com/moler/index_ncm.html}. \cr
//' Kai Habel, Jaroslav Hajek, interp2, Octave.
//' @examples
//' # example in MatLab
//' library(lattice)
//' # data generation
//' x <- seq(-3, 3, 1)
//' xm <- expand.grid(x, x)
//' z <- 3*(1-xm[,1])^2.*exp(-(xm[,1]^2) - (xm[,2]+1)^2) -
//'   10*(xm[,1]/5 - xm[,1]^3 - xm[,2]^5)*exp(-xm[,1]^2-xm[,2]^2) -
//'   1/3*exp(-(xm[,1]+1)^2 - xm[,2]^2)
//' dat <- data.frame(xm, z)
//' # graph of original data
//' wireframe(z ~ Var1 + Var2, dat, drape = TRUE, colorkey = TRUE)
//'
//' xi <- seq(-3, 3, 0.25)
//' zi_l <- interp2(x, x, matrix(z, length(x)), xi, xi, 'linear')
//' dat_l <- data.frame(expand.grid(xi, xi), zi_l)
//' # graph of linearly interpolation
//' wireframe(zi_l ~ Var1 + Var2, dat_l, drape = TRUE, colorkey = TRUE)
//'
//' zi_s <- interp2(x, x, matrix(z, length(x)), xi, xi, 'spline')
//' dat_s <- data.frame(expand.grid(xi, xi), zi_s)
//' # graph of interpolation with spline
//' wireframe(zi_s ~ Var1 + Var2, dat_s, drape = TRUE, colorkey = TRUE)
//' @export
// [[Rcpp::export]]
arma::mat interp2(const arma::vec& x, const arma::vec& y, const arma::mat& v,
                  const arma::vec& xi, arma::vec& yi, const std::string& method){
  chk_mat(x, "x");
  chk_mat(y, "y");
  chk_mat(v, "v");
  chk_mat(xi, "xi");
  chk_mat(yi, "yi");
  if (x.n_elem != v.n_cols)
    Rcpp::stop("The number of columns of v must be equal to the length of x.");
  if (y.n_elem != v.n_rows)
    Rcpp::stop("The number of rows of v must be equal to the length of y.");
  if (!x.is_sorted())
    RMessage("x are not strictly monotonic increasing.\nx will be sorted.");
  if (method != "linear" && method != "spline")
    Rcpp::stop("Method only support linear and spline.\n");

  uvec xu_idx = find_unique(x);
  vec xu = x(xu_idx);
  mat v_tmp = v.cols(xu_idx);
  if (!xu.is_sorted()) {
    uvec si = sort_index(xu);
    xu = xu(si);
    v_tmp = v_tmp.cols(si);
  }
  if (xu.n_elem != x.n_elem) {
    RMessage("The grid vectors are not strictly monotonic increasing.");
    RMessage("The values of v for duplicated values of x will be averaged.");
    for (uword k = 0; k < xu.n_elem; ++k)
      v_tmp.col(k) = mean(v.cols(find(x == xu(k))), 1);
  }

  if (!y.is_sorted())
    RMessage("y are not strictly monotonic increasing.\ny will be sorted.");
  uvec yu_idx = find_unique(y);
  vec yu = y(yu_idx);
  mat vu = v_tmp.rows(yu_idx);
  if (!yu.is_sorted()) {
    uvec si = sort_index(yu);
    yu = yu(si);
    vu = vu.rows(si);
  }
  if (yu.n_elem != y.n_elem) {
    RMessage("The grid vectors are not strictly monotonic increasing.");
    RMessage("The values of v for duplicated values of y will be averaged.");
    for (uword k = 0; k < yu.n_elem; ++k)
      vu.row(k) = mean(v_tmp.rows(find(y == yu(k))));
  }

  mat vi(xi.n_elem, yi.n_elem);
  if (method == "linear") {
    uvec xidx_tmp = lookup(xu, xi), yidx_tmp = lookup(yu, yi);
    xidx_tmp(find(xidx_tmp == xu.n_elem)).fill(xu.n_elem - 1);
    yidx_tmp(find(yidx_tmp == yu.n_elem)).fill(yu.n_elem - 1);
    xidx_tmp--;
    yidx_tmp--;

    umat xidx = repmat(xidx_tmp.t(), yi.n_elem, 1),
      yidx = repmat(yidx_tmp, 1, xi.n_elem);

    uword nvr = vu.n_rows, nvc = vu.n_cols;
    mat a = vu.submat(0, 0, nvr-2, nvc-2),
      b = vu.submat(0, 1, nvr-2, nvc-1) - a,
      c = vu.submat(1, 0, nvr-1, nvc-2) - a,
      d = vu.submat(1, 1, nvr-1, nvc-1) - a - b - c;

    vec dx = diff(xu), dy = diff(yu);
    vec xsc = (vectorise(repmat(xi.t(), yi.n_elem, 1)) - xu.elem(vectorise(xidx))) / dx.elem(vectorise(xidx)),
      ysc = (repmat(yi, xi.n_elem, 1) - yu.elem(vectorise(yidx))) / dy.elem(vectorise(yidx));
    uvec idx = vectorise(yidx) + a.n_rows*vectorise(xidx);
    vi = reshape(a(idx) + b(idx) % xsc + c(idx) % ysc + d(idx) % xsc % ysc,
                 yi.n_elem, xi.n_elem);
  } else if (method == "spline") {
    vi = spline_cpp(yu, vu, yi);
    vi = spline_cpp(xu, vi.t(), xi).t();
  }
  return vi;
}
