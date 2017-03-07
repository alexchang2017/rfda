#ifndef getFpcScoresReg_h
#define getFpcScoresReg_h

#include "common.h"
Rcpp::List getFpcScoresReg(const arma::vec& allTimePnts, const arma::vec& splitVar, const arma::mat& yMat,
                           const arma::mat& eigFuncs, const arma::rowvec& eigVals, const arma::vec& measErrVar,
                           const std::string& methodFPCS);

#endif
