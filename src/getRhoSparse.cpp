#include "common.h"
#include "getFpcScoresSparse.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

// [[Rcpp::export]]
arma::vec getRhoSparse(const arma::vec& rhoFactor, const arma::vec& splitVar, const arma::field<arma::vec>& yList,
                       const arma::vec& timeIdx, const arma::vec& subId, const arma::mat& eigFuncs, const arma::rowvec& eigVals,
                       const arma::vec& measErrVar, const std::string& methodFPCS, const std::string& rho){
  vec measErrVar2 = measErrVar;
  Rcpp::List fpcScoresResTmp;
  for (uword i = 0; i < 2; ++i) {
    fpcScoresResTmp = getFpcScoresSparse(splitVar, yList, timeIdx, subId, eigFuncs, eigVals, measErrVar2, methodFPCS, 1.0);
    measErrVar2 = mean(Rcpp::as<mat>(fpcScoresResTmp[2])).t();
  }

  // uvec cvIdx = zeros<uvec>(yList.n_elem);
  // if (rho == "cv-random") {
  //   for (uword i = 0; i < yList.n_elem; ++i)
  //     cvIdx(i) = as_scalar(randi<uvec>(1, distr_param(0, yList(i).n_elem - 1)));
  // } else {
  //   for (uword i = 0; i < yList.n_elem; ++i)
  //     cvIdx(i) = (1001 + i) % yList(i).n_elem;
  // }

  // vec rhoCand = rhoFactor * linspace<vec>(0.01, 0.23, 51);
  // rhoCand.elem(find(rhoCand > ))

  // optnsTmp <- optns
  // optnsTmp$verbose <- FALSE
  // for (j in 1:2) {
  //   yhat <- GetCEScores(y, t, optnsTmp, mu, obsGrid, fittedCov, lambda, phi, sigma2)[3, ]
  //   sigma2 <- mean(mapply(function(a, b) mean((a - b)^2, na.rm=TRUE), yhat, y), na.rm=TRUE)
  // }
  //
  // R <- sqrt((trapzRcpp(obsGrid, mu ^ 2) + sum(lambda)) / diff(range(obsGrid)))
  //   a1 <- 0.01; a2 <- 0.22
  //   etaCand <- seq(a1, a2, length.out=50)
  //     rhoCand <- etaCand * R
  //     rhoCand <- rhoCand[rhoCand > sigma2]
  //   rhoCand <- c(sigma2, rhoCand)
  //
  //     leaveOutInd <- RandTime(t, isRandom=FALSE)
  //
  //     cvScores <- sapply(rhoCand, cvRho, leaveOutInd=leaveOutInd, y=y, t=t, optns=optns, mu=mu, obsGrid=obsGrid, fittedCov=fittedCov, lambda=lambda, phi=phi)
  //
  //

  // vec measErrVarVec = measErrVar;
  // mat fpcScores = zeros<mat>(eigFuncs.n_cols, subId.n_elem);
  // field<mat> fpcScoresVar(uniSubId.n_elem);

  // WorkerFpcScoresSparse FpcScoresSparse_worker(uniVars.n_elem, uniSubId, subId2, timeIdx2, measErrVarVec, eigVals, yList,
  //                                              eigFuncs, methodFPCS, fpcScores, fpcScoresVar);
  // RcppParallel::parallelFor(0, uniSubId.n_elem, FpcScoresSparse_worker);
  // inplace_trans(fpcScores);

  return measErrVar2;
}

