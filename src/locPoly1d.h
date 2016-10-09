#ifndef locPoly1d_H
#define locPoly1d_H

#include "common.h"
Rcpp::NumericVector locPoly1d(const double& bandwidth, const arma::vec& x, const arma::vec& y,
                              const arma::vec& w, const arma::vec& xout, const std::string& kernel,
                              const double& drv, const double& degree);
double gcvLocPoly1d(arma::vec bwCand, const arma::vec& x, const arma::vec& y,
                    const arma::vec& w, const std::string& kernel,
                    const double& drv, const double& degree);

#endif
