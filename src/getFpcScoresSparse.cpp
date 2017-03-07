#include "common.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct WorkerFpcScoresSparse: public RcppParallel::Worker {
  const uword& numVar;
  const uvec& uniSubId;
  const uvec& subId2;
  const uvec& timeIdx2;
  const vec& measErrVarVec;
  const rowvec& eigVals;
  const field<vec>& yList;
  const mat& eigFuncs;
  const std::string& methodFPCS;
  const double& getMse;
  mat& fpcScores;
  field<mat>& fpcScoresVar;
  mat& mse;

  WorkerFpcScoresSparse(const uword& numVar, const uvec& uniSubId, const uvec& subId2, const uvec& timeIdx2,
                        const vec& measErrVarVec, const rowvec& eigVals, const field<vec>& yList, const mat& eigFuncs,
                        const std::string& methodFPCS, const double& getMse, mat& fpcScores, field<mat>& fpcScoresVar, mat& mse):
    numVar(numVar), uniSubId(uniSubId), subId2(subId2), timeIdx2(timeIdx2), measErrVarVec(measErrVarVec), eigVals(eigVals),
    yList(yList), eigFuncs(eigFuncs), methodFPCS(methodFPCS), getMse(getMse), fpcScores(fpcScores), fpcScoresVar(fpcScoresVar),
    mse(mse) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (uword i = begin; i < end; ++i) {
      uvec idxSub = find(subId2 == uniSubId(i));
      uvec idxTpt = timeIdx2(idxSub);
      uvec mappingIdx = repmat(idxTpt, numVar, 1);
      mat eigFuncsSub = eigFuncs.rows(mappingIdx);
      vec measErrVarVecSub = measErrVarVec.elem(mappingIdx);

      if (methodFPCS == "CE") {
        mat tmp = eigFuncsSub.each_row() % eigVals;
        inplace_trans(tmp);
        mat fiitedCovWithMeVar = eigFuncsSub * tmp;
        fiitedCovWithMeVar.diag() += measErrVarVecSub;
        mat fiitedCovWithMeVarInv = pinv(fiitedCovWithMeVar);
        fpcScores.col(i) = tmp * fiitedCovWithMeVarInv * yList(i);
        fpcScoresVar(i) = -tmp * fiitedCovWithMeVarInv * tmp.t();
        if (getMse) {
          for (uword j = 0; j < numVar; ++j)
            mse(i, j) = mean(square(yList(i).subvec(idxTpt.n_elem * j, idxTpt.n_elem * (j+1) - 1) -
              eigFuncsSub.rows(idxTpt.n_elem * j, idxTpt.n_elem * (j+1) - 1) * fpcScores.col(i)));
        }
      } else if (methodFPCS == "LS") {
        mat Ai = pinv(eigFuncsSub.t() * eigFuncsSub);
        mat tmp = Ai * eigFuncsSub.t();
        fpcScores.col(i) = tmp * yList(i);
        fpcScoresVar(i) = (tmp.each_row() % measErrVarVecSub.t()) * tmp.t();
        if (getMse) {
          for (uword j = 0; j < numVar; ++j)
            mse(i, j) = mean(square(yList(i).subvec(idxTpt.n_elem * j, idxTpt.n_elem * (j + 1) - 1) -
              eigFuncsSub.rows(idxTpt.n_elem * j, idxTpt.n_elem * (j + 1) - 1) * fpcScores.col(i)));
        }
      } else if (methodFPCS == "WLS") {
        mat tmp = eigFuncsSub.each_col() % pow(measErrVarVecSub, -1);
        inplace_trans(tmp);
        mat tmp2 = pinv(tmp * eigFuncsSub);
        fpcScores.col(i) = tmp2 * tmp * yList(i);
        fpcScoresVar(i) = tmp2;
        if (getMse) {
          for (uword j = 0; j < numVar; ++j)
            mse(i, j) = mean(square(yList(i).subvec(idxTpt.n_elem * j, idxTpt.n_elem * (j+1) - 1) -
              eigFuncsSub.rows(idxTpt.n_elem * j, idxTpt.n_elem * ( j+ 1) - 1) * fpcScores.col(i)));
        }
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::List getFpcScoresSparse(const arma::vec& splitVar, const arma::field<arma::vec>& yList, const arma::vec& timeIdx,
                              const arma::vec& subId, const arma::mat& eigFuncs, const arma::rowvec& eigVals,
                              const arma::vec& measErrVar, const std::string& methodFPCS, const double& getMse){
  chk_mat(splitVar, "splitVar");
  chk_mat(timeIdx, "timeIdx");
  chk_mat(subId, "subId");
  chk_mat(eigFuncs, "eigFuncs");
  chk_mat(eigVals, "eigVals");
  chk_mat(measErrVar, "measErrVar");

  uvec subId2 = conv_to<uvec>::from(subId - 1);
  uvec timeIdx2 = conv_to<uvec>::from(timeIdx - 1);
  uvec idx = conv_to<uvec>::from(splitVar - 1.0);
  uvec uniVars = sort(unique(idx));
  uvec uniSubId = sort(unique(subId2));
  if (uniSubId.n_elem != yList.n_rows)
    Rcpp::stop("The length of yList must equal to the number of unique subId!");
  if (timeIdx.n_elem != subId.n_elem)
    Rcpp::stop("The length of timeIdx must equal to the length of subId!");
  if (splitVar.n_elem != eigFuncs.n_rows)
    Rcpp::stop("The length of splitVar must be equal to the number of rows of eigFuncs and yMat!");
  if (uniVars.n_elem != measErrVar.n_elem)
    Rcpp::stop("The length of measErrVar must be equal to the number of unique splitVar!");
  if (any(measErrVar <= 0))
    Rcpp::stop("The elements of measErrVar must be positive!");
  if (eigFuncs.n_cols != eigVals.n_elem)
    Rcpp::stop("The number of columns of eigFuncs must be equal to the length of eigVals!");
  if (!is_finite(getMse) || (getMse != 0.0 && getMse != 1.0))
    Rcpp::stop("mse must be 0 or 1.\n");

  vec measErrVarVec = measErrVar.elem(idx);
  mat mse = zeros<mat>(uniSubId.n_elem, uniVars.n_elem);
  mat fpcScores = zeros<mat>(eigFuncs.n_cols, uniSubId.n_elem);
  field<mat> fpcScoresVar(uniSubId.n_elem);
  for (uword i = 0; i < uniSubId.n_elem; ++i)
    fpcScoresVar(i) = zeros<mat>(eigFuncs.n_cols, eigFuncs.n_cols);

  WorkerFpcScoresSparse FpcScoresSparse_worker(uniVars.n_elem, uniSubId, subId2, timeIdx2, measErrVarVec, eigVals, yList,
                                               eigFuncs, methodFPCS, getMse, fpcScores, fpcScoresVar, mse);
  RcppParallel::parallelFor(0, uniSubId.n_elem, FpcScoresSparse_worker);
  inplace_trans(fpcScores);

  if (getMse)
    return Rcpp::List::create(Rcpp::Named("fpcScores") = fpcScores,
                              Rcpp::Named("fpcScoresVar") = fpcScoresVar,
                              Rcpp::Named("fpcScoresVar") = mse);

  return Rcpp::List::create(Rcpp::Named("fpcScores") = fpcScores,
                            Rcpp::Named("fpcScoresVar") = fpcScoresVar,
                            Rcpp::Named("fpcScoresVar") = R_NilValue);
}

