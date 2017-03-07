#include "common.h"

// [[Rcpp::export]]
Rcpp::List getFpcScoresReg(const arma::vec& allTimePnts, const arma::vec& splitVar, const arma::mat& yMat,
                           const arma::mat& eigFuncs, const arma::rowvec& eigVals, const arma::vec& measErrVar,
                           const std::string& methodFPCS){
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

  uvec idx = conv_to<uvec>::from(splitVar - 1.0);
  vec measErrVarVec = measErrVar.elem(idx);
  mat fpcScores = zeros<mat>(eigFuncs.n_cols, yMat.n_cols);
  mat fpcScoresVar = zeros<mat>(eigFuncs.n_cols, eigFuncs.n_cols);
  if (methodFPCS == "CE") {
    mat tmp = eigFuncs.each_row() % eigVals;
    inplace_trans(tmp);
    mat fiitedCovWithMeVar = eigFuncs * tmp;
    fiitedCovWithMeVar.diag() += measErrVarVec;
    mat fiitedCovWithMeVarInv = pinv(fiitedCovWithMeVar);
    fpcScores = tmp * fiitedCovWithMeVarInv * yMat;
    inplace_trans(fpcScores);
    fpcScoresVar = -tmp * fiitedCovWithMeVarInv * tmp.t();
  } else if (methodFPCS == "LS") {
    mat Ai = pinv(eigFuncs.t() * eigFuncs);
    mat tmp = Ai * eigFuncs.t();
    fpcScores = tmp * yMat;
    inplace_trans(fpcScores);
    fpcScoresVar = (tmp.each_row() % measErrVarVec.t()) * tmp.t();
  } else if (methodFPCS == "WLS") {
    mat tmp = eigFuncs.each_col() % pow(measErrVarVec, -1);
    inplace_trans(tmp);
    mat tmp2 = pinv(tmp * eigFuncs);
    fpcScores = tmp2 * tmp * yMat;
    fpcScoresVar = tmp2;
  }
  return Rcpp::List::create(Rcpp::Named("fpcScores") = fpcScores,
                            Rcpp::Named("fpcScoresVar") = fpcScoresVar);
}

