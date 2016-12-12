#include "common.h"
#include "interp.h"
#include <string>

// [[Rcpp::export]]
Rcpp::List getEigRes(const arma::mat& CFMat2, const arma::vec& variable, const arma::vec& workTimePnts,
                     const arma::vec& meanFuncsWork, const arma::vec& allTimePnts){
  chk_mat(CFMat2, "CFMat2", "double");
  chk_mat(variable, "variable", "double");
  chk_mat(workTimePnts, "workTimePnts", "double");
  chk_mat(meanFuncsWork, "meanFuncsWork", "double");
  chk_mat(allTimePnts, "allTimePnts", "double");

  vec uniVars = sort(unique(variable));
  if (!CFMat2.is_square())
    Rcpp::stop("CFMat2 must be square matrix!");
  if (CFMat2.n_rows != variable.n_elem || meanFuncsWork.n_elem != variable.n_elem)
    Rcpp::stop("The lengths of variable and meanFuncsWork must be eqaul to the number of rows of CFMat2!");
  if (uniVars.n_elem * workTimePnts.n_elem != meanFuncsWork.n_elem)
    Rcpp::stop("The length of meanFuncsWork must equal to the length of workTimePnts multiply the number of unique variable!");

  // Perform eigen decomposition
  vec eigVals;
  mat eigFuncsWork;
  eig_sym(eigVals, eigFuncsWork, CFMat2);

  // remove nonpositive eigenvalues and the corresponding eigenfunctions
  uvec idx = find(eigVals > 0);
  uword numEigs = CFMat2.n_cols - idx.n_elem;
  Rcpp::Function RMessage2("message");
  if (any(eigVals <= 0.0))
    RMessage2("Warning: ", numEigs, " real eigenvalues are negative or zero and are removed!");
  eigVals = eigVals.elem(idx);
  eigFuncsWork = eigFuncsWork.cols(idx);

  // sort the eigenvalues and the corresponding eigenfunctions
  double tpAvgDiff = (workTimePnts.max() - workTimePnts.min()) / ((double) workTimePnts.n_elem - 1.0);
  uvec sortIdx = sort_index(eigVals, "descend");
  eigVals = eigVals.elem(sortIdx) * tpAvgDiff;
  eigFuncsWork = eigFuncsWork.cols(sortIdx);

  // normalize the eigenfunctions
  rowvec normValsWork = zeros<rowvec>(eigFuncsWork.n_cols);
  for (uword i = 0; i < uniVars.n_elem; i++)
    normValsWork += trapz_cpp(workTimePnts, square(eigFuncsWork.rows(find(variable == uniVars(i)))));
  eigFuncsWork.each_row() /= sqrt(normValsWork);
  eigFuncsWork.each_row() %= sign(cov(meanFuncsWork, eigFuncsWork));

  // get eigenfunctions onto allTimePnts
  std::string interp1_method = "spline";
  mat eigFuncs = zeros<mat>(allTimePnts.n_elem * uniVars.n_elem, eigFuncsWork.n_cols);
  rowvec normVals = zeros<rowvec>(eigFuncsWork.n_cols);
  for (uword i = 0; i < uniVars.n_elem; i++) {
    eigFuncs.rows(i * allTimePnts.n_elem , (i + 1) * allTimePnts.n_elem - 1) =
      interp1_cpp(workTimePnts, eigFuncsWork.rows(find(variable == uniVars(i))), allTimePnts, interp1_method);
    normVals += trapz_cpp(allTimePnts, square(eigFuncs.rows(i * allTimePnts.n_elem , (i + 1) * allTimePnts.n_elem - 1)));
  }
  eigFuncs.each_row() /= sqrt(normVals);

  Rcpp::NumericVector eigValsR = Rcpp::wrap(eigVals);
  eigValsR.attr("dim") = R_NilValue;
  return Rcpp::List::create(Rcpp::Named("eigVals") = eigValsR,
                            Rcpp::Named("eigFuncsWork") = eigFuncsWork,
                            Rcpp::Named("eigFuncs") = eigFuncs);
}
