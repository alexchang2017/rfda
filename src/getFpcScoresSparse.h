#ifndef getFpcScoresSparse_h
#define getFpcScoresSparse_h

#include "common.h"
Rcpp::List getFpcScoresSparse(const arma::vec& splitVar, const arma::field<arma::vec>& yList, const arma::vec& timeIdx,
                              const arma::vec& subId, const arma::mat& eigFuncs, const arma::rowvec& eigVals,
                              const arma::vec& measErrVar, const std::string& methodFPCS, const double& getMse);

#endif
