#ifndef interp_H
#define interp_H

#include "common.h"
arma::mat spline_cpp(const arma::vec& x, const arma::mat& y, const arma::vec& xi);
arma::mat interp1_cpp(const arma::vec& x, const arma::mat& y, const arma::vec& xi,
                  const std::string& method);
arma::mat interp2(const arma::vec& x, const arma::vec& y, const arma::mat& v,
                  const arma::vec& xi, arma::vec& yi, const std::string& method);

#endif
