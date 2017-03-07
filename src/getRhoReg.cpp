#include "common.h"
#include "getFpcScoresReg.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

// struct cvRhoReg_Worker: public Worker {
//   const mat& phiSub;
//   const vec& muSub;
//   const rowvec& lambdaWeight;
//   const vec& rho_candidates;
//   const double& sigma1;
//   const uvec& tjID;
//   const mat& lambdaMat;
//   const field<vec>& yField;
//   const std::string& estFPCscores;
//   const uvec& ni;
//   vec& cv;
//   cvrhoWorker_regular2(const mat& phiSub, const vec& muSub, const rowvec& lambdaWeight,
//                        const vec& rho_candidates, const double& sigma1, const uvec& tjID,
//                        const mat& lambdaMat, const field<vec>& yField,
//                        const std::string& estFPCscores, const uvec& ni, vec& cv):
//     phiSub(phiSub), muSub(muSub), lambdaWeight(lambdaWeight), rho_candidates(rho_candidates),
//     sigma1(sigma1), tjID(tjID), lambdaMat(lambdaMat), yField(yField), estFPCscores(estFPCscores),
//     ni(ni), cv(cv) {}
//   void operator()(std::size_t begin, std::size_t end)
//   {
//     mat Atmp, D, A, Ai;
//     rowvec xi_est;
//     uvec sampleID;
//     double sigmaTmp = sigma1;
//     for (uword k = begin; k < end; ++k)
//     {
//       sigmaTmp = (sigmaTmp >= rho_candidates(k)) ? sigma1 : rho_candidates(k);
//       for (uword i = 0; i < yField.n_elem; ++i)
//       {
//         sampleID = linspace<uvec>(0, ni(i)-1, ni(i));
//         sampleID = sampleID.elem(find(sampleID != tjID(i)));
//         if (estFPCscores.compare("CE") == 0)
//         {
//           Atmp = lambdaMat * phiSub.rows(sampleID).t();
//           D = phiSub.rows(sampleID) * Atmp;
//           D.diag() += sigmaTmp;
//           A = Atmp * pinv(D);
//           xi_est = (A * (yField(i).elem(sampleID) - muSub.elem(sampleID))).t();
//         } else if (estFPCscores.compare("LS") == 0)
//         {
//           A = phiSub.rows(sampleID).t() * phiSub.rows(sampleID);
//           A = (A + A.t()) / 2;
//           Ai = pinv(A);
//           xi_est = (yField(i).elem(sampleID) - muSub.elem(sampleID)).t() * phiSub.rows(sampleID) * Ai.t();
//         } else if (estFPCscores.compare("WLS") == 0)
//         {
//           Atmp = sigmaTmp * eye<mat>(ni(0)-1, ni(0)-1);
//           Ai = pinv(Atmp);
//           A = phiSub.rows(sampleID).t() * Ai * phiSub.rows(sampleID);
//           A = (A + A.t()) / 2;
//           xi_est = (yField(i).elem(sampleID) - muSub.elem(sampleID)).t() * Ai.t() * phiSub.rows(sampleID) * pinv(A).t();
//         }
//         cv(k) += std::pow(muSub(tjID(i)) + dot(xi_est, phiSub.row(tjID(i))) - yField(i)(tjID(i)), 2.0);
//       }
//     }
//   }
// };

// struct Worker_getLogLikMilr : public RcppParallel::Worker {
//   const uvec& bag2;
//   const vec& y;
//   const mat& X;
//   const vec& beta;
//   double logLikMilr;
//
//   Worker_getLogLikMilr(const uvec& bag2, const vec& y, const mat& X, const vec& beta):
//     bag2(bag2), y(y), X(X), beta(beta), logLikMilr(0) {}
//
//   Worker_getLogLikMilr(const Worker_getLogLikMilr& getLogLikMilr_worker, RcppParallel::Split):
//     bag2(getLogLikMilr_worker.bag2), y(getLogLikMilr_worker.y), X(getLogLikMilr_worker.X),
//     beta(getLogLikMilr_worker.beta), logLikMilr(0) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for (uword i = begin; i < end; ++i) {
//       uvec idx = find(bag2 == i);
//       double prob = 1 - prod(1.0 - logit(X.rows(idx), beta));
//       logLikMilr += y(idx(0)) * log(prob) + (1 - y(idx(0))) * log(1 - prob);
//     }
//   }
//
//   // join my value with that of another Sum
//   void join(const Worker_getLogLikMilr& rhs) {
//     logLikMilr += rhs.logLikMilr;
//   }
// };

// [[Rcpp::export]]
arma::vec getRhoReg(const arma::vec& rhoFactor, const arma::vec& allTimePnts, const arma::vec& splitVar,
                    const arma::mat& yMat, const arma::mat& eigFuncs, const arma::rowvec& eigVals,
                    const arma::vec& measErrVar, const std::string& methodFPCS, const std::string& rho,
                    const double& minMearErr) {
  vec measErrVar2 = measErrVar;
  vec tmp = zeros<vec>(splitVar.n_elem);
  vec uniVars = sort(unique(splitVar));
  Rcpp::List fpcScoresResTmp;
  for (uword i = 0; i < 2; ++i) {
    fpcScoresResTmp = getFpcScoresReg(allTimePnts, splitVar, yMat, eigFuncs, eigVals, measErrVar2, methodFPCS);
    tmp = mean(square(eigFuncs * Rcpp::as<mat>(fpcScoresResTmp[0]) - yMat), 1);
    for (uword j = 0; j < uniVars.n_elem; ++j)
      measErrVar2(j) = mean(tmp.elem(find(splitVar == uniVars(j))));
    measErrVar2.elem(find(measErrVar2 < minMearErr)).fill(minMearErr);
  }

  umat testIdx = zeros<umat>(yMat.n_cols, uniVars.n_elem);
  if (rho == "cv-random") {
    testIdx.col(0) = randi<uvec>(yMat.n_cols, distr_param(0, allTimePnts.n_elem - 1));
  } else {
    uvec tmpIdx = linspace<uvec>(1000, 1000 + yMat.n_cols - 1, yMat.n_cols);
    testIdx.col(0) = tmpIdx - floor(tmpIdx / allTimePnts.n_elem) * allTimePnts.n_elem;
  }

  if (uniVars.n_elem > 1) {
    for (uword j = 1; j < uniVars.n_elem; ++j)
      testIdx.col(j) = testIdx.col(j - 1) + allTimePnts.n_elem;
  }

  vec cvScores = zeros<vec>(uniVars.n_elem);
  uvec allIdx = linspace<uvec>(0, yMat.n_rows - 1, yMat.n_rows);

  vec grid = linspace<vec>(0.01, 0.22, 50);
  field<vec> rhoCandList(uniVars.n_elem);
  for (uword j = 0; j < uniVars.n_elem; ++j)
    rhoCandList(j) = join_cols(measErrVar2(j) * ones<vec>(1),
                rhoFactor(j) * grid.elem(find(grid * rhoFactor(j) > measErrVar2(j))));

  // WorkerCvRho cvRho_worker();
  // RcppParallel::parallelFor(0, yMat.n_cols, cvRho_worker);

  return measErrVar2;
}
