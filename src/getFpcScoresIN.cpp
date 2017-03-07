#include "common.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct WorkerFpcScoresIN: public RcppParallel::Worker {
  const vec& allTimePnts;
  const umat& idxMat;
  const mat& yMat;
  const mat& eigFuncs;
  const mat& shrinkFactor;
  mat& fpcScores;

  WorkerFpcScoresIN(const vec& allTimePnts, const umat& idxMat, const mat& yMat, const mat& eigFuncs,
                    const mat& shrinkFactor, mat& fpcScores):
    allTimePnts(allTimePnts), idxMat(idxMat), yMat(yMat), eigFuncs(eigFuncs),
    shrinkFactor(shrinkFactor), fpcScores(fpcScores) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (uword i = begin; i < end; ++i) {
      for (uword j = 0; j < idxMat.n_cols; ++j) {
        mat tmpEigFuncs = eigFuncs.rows(idxMat.col(j));
        fpcScores.row(i) += trapz_cpp(allTimePnts, tmpEigFuncs.each_col() % yMat(idxMat.col(j), i * ones<uvec>(1))) %
          shrinkFactor.row(j);
      }
    }
  }
};

// [[Rcpp::export]]
arma::mat getFpcScoresIN(const arma::vec& allTimePnts, const arma::vec& splitVar, const arma::mat& yMat,
                         const arma::mat& eigFuncs, const double& shrink, const arma::rowvec& eigVals,
                         const arma::vec& measErrVar){
  chk_mat(allTimePnts, "allTimePnts");
  chk_mat(splitVar, "splitVar");
  chk_mat(yMat, "yMat");
  chk_mat(eigFuncs, "eigFuncs");
  chk_mat(eigVals, "eigVals");
  chk_mat(measErrVar, "measErrVar");

  vec uniVars = sort(unique(splitVar));
  if (allTimePnts.n_elem * uniVars.n_elem != eigFuncs.n_rows)
    Rcpp::stop("The number of rows of eigFuncs must equal to the length of allTimePnts multiply the number of unique splitVar!");
  if (splitVar.n_elem != eigFuncs.n_rows || splitVar.n_elem != yMat.n_rows)
    Rcpp::stop("The length of splitVar must be equal to the number of rows of eigFuncs and yMat!");
  if (uniVars.n_elem != measErrVar.n_elem)
    Rcpp::stop("The length of measErrVar must be equal to the number of unique splitVar!");
  if (any(measErrVar <= 0))
    Rcpp::stop("The elements of measErrVar must be positive!");
  if (eigFuncs.n_cols != eigVals.n_elem)
    Rcpp::stop("The number of columns of eigFuncs must be equal to the length of eigVals!");
  if (std::abs(shrink - std::floor(shrink)) > 1e-6 || (shrink != 0.0 && shrink != 1.0))
    Rcpp::stop("shrink must be TRUE or FALSE.");

  mat shrinkFactor = ones<mat>(measErrVar.n_elem, eigVals.n_elem);
  if (shrink == 1.0) {
    double timeRange = allTimePnts.max() - allTimePnts.min();
    for (uword j = 0; j < uniVars.n_elem; ++j)
      shrinkFactor.row(j) = eigVals / (eigVals + timeRange * measErrVar(j) / (double) allTimePnts.n_elem);
  }

  umat idxMat = zeros<umat>(allTimePnts.n_elem, uniVars.n_elem);
  for (uword j = 0; j < uniVars.n_elem; ++j)
    idxMat.col(j) = find(splitVar == uniVars(j));

  mat fpcScores = zeros<mat>(yMat.n_cols, eigFuncs.n_cols);
  WorkerFpcScoresIN FpcScoresIN_worker(allTimePnts, idxMat, yMat, eigFuncs, shrinkFactor, fpcScores);
  RcppParallel::parallelFor(0, yMat.n_cols, FpcScoresIN_worker);
  return fpcScores;
}

