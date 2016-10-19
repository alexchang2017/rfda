#include "common.h"
#include <cmath>

//' Correlation Generation Function
//'
//' This function performs the generation of correlation sequence.
//'
//' The Matern correlation is \deqn{rho1 = 2 * sqrt(n) * rho / rho0} and
//' \deqn{r = rho1 ^ nu * besselK(nu, rho1) / gamma(nu) / 2 ^ (nu - 1)}
//' As nu goes to infinity, the correlation converges to \deqn{exp(-(rho/rho0)^2)}.
//'
//' @param x A numerical vector.
//' @param corrType A string, the type of correlation function, "BesselJ", "Matern" or "rq".
//'   Please see "Details" for more details.
//' @param x0 A numeric, the shape of correlation. It must be greater than zero. Default value is 1.
//' @param nu A numeric, it is only used when corrType is 'Matern'. It must be greater than zero.
//'   Default value is 2.5. Please see "Details" for more details.
//' @return A correlation sequence.
//' @examples
//' require(ggplot2)
//' x <- seq(-10, 10, len = 101)
//' ggplot(data.frame(x, y = corrGen(x, "BesselJ")), aes(x, y, colour = "BesselJ")) +
//'   geom_line(size = 1.2) + geom_line(aes(x, y, colour = "Matern (nu=2.5)"),
//'   data.frame(x, y = corrGen(x, "Matern")), size = 1.2) +
//'   geom_line(aes(x, y, colour = "rq"), data.frame(x, y = corrGen(x, "rq")), size = 1.2) +
//'   labs(title = "Three types of correlation sequences", colour = "Correlation Type",
//'        x = "Distance", y = "Correlation Strength")
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector corrGen(const arma::vec& x, const std::string& corrType,
                            const double x0 = 1.0, const double nu = 2.5){
  // data checking
  chk_mat(x, "x", "double");
  if (!is_finite(x0) || x0 <= 0)
    Rcpp::stop("x0 must be a positive number.\n");
  if (!is_finite(nu) || nu <= 0)
    Rcpp::stop("nu must be a positive number.\n");
  if (corrType != "BesselJ" && corrType != "Matern" && corrType != "rq")
    Rcpp::stop("Unsupported corrType. corrType must be 'BesselJ', 'Matern' or 'rq'.\n");

  // import R function
  Rcpp::Function RbesselJ("besselJ");
  Rcpp::Function RbesselK("besselK");
  Rcpp::Function Rgamma("gamma");

  // initialize data
  Rcpp::NumericVector outVec(x.n_elem);
  if (corrType == "BesselJ"){
    outVec = RbesselJ(Rcpp::wrap(abs(x) / x0), 0.0);
  } else if (corrType == "Matern") {
    double c = Rcpp::as<double>(Rgamma(nu)) * std::pow(2.0, nu-1.0);
    vec x1 = 2.0 * std::sqrt(nu) * abs(x) / x0;
    vec bk = Rcpp::as<vec>(RbesselK(Rcpp::wrap(x1), nu));
    vec out = pow(x1, nu) % bk / c;
    out.elem(find_nonfinite(out)).ones();
    outVec = Rcpp::wrap(out);
  } else if (corrType == "rq") {
    outVec = Rcpp::wrap(pow(1.0 + square(x / x0), -1.0));
  }
  // remove dimesion of output
  outVec.attr("dim") = R_NilValue;
  return(outVec);
}

//' Correlation Generation Function
//'
//' This function performs the generation of correlation sequence.
//'
//' The Matern correlation is \deqn{rho1 = 2 * sqrt(n) * rho / rho0} and
//' \deqn{r = rho1 ^ nu * besselK(nu, rho1) / gamma(nu) / 2 ^ (nu - 1)}
//' As nu goes to infinity, the correlation converges to \deqn{exp(-(rho/rho0)^2)}.
//'
//' @param n An integer, the number of sample size.
//' @param timePnt A numeric vector, the observed time.
//' @param meanFunc,varFunc The function to generate mean function and variance function.
//' @param measErrVar The variance of measurement error.
//' @param corrType,x0,nu The parameters to generate correlation sequence.
//' @return A data.frame containing sample id, observed time and corresponding observed values.
//' @examples
//' require(ggplot2)
//' tp <- seq(1, 10, len = 7)
//' DT <- funcDataGen(6, tp, function(x) sin(x), function(x) rep(1, length(x)), "BesselJ")
//' ggplot(DT, aes(x = t, y = y, color = factor(sampleID))) + geom_line() + labs(color = "Sample ID")
//'
//' DT2 <- funcDataGen(6, tp, function(x) sin(x), function(x) cos(x)/2 + 2, "BesselJ")
//' ggplot(DT2, aes(x = t, y = y, color = factor(sampleID))) + geom_line() + labs(color = "Sample ID")
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame funcDataGen(const double& n, const arma::vec& timePnt, const Rcpp::Function meanFunc,
                            const Rcpp::Function varFunc, const std::string corrType,
                            const double measErrVar = 1, const double x0 = 1.0, const double nu = 2.5){
  // data checking
  chk_mat(timePnt, "timePnt", "double");
  if (!is_finite(n) || n <= 0 || std::abs(n - std::floor(n)) > 1e-6)
    Rcpp::stop("n must be a positive integer.\n");
  if (corrType != "BesselJ" && corrType != "Matern" && corrType != "rq")
    Rcpp::stop("Unsupported corrType. corrType must be 'BesselJ', 'Matern' or 'rq'.\n");
  if (!is_finite(measErrVar) || measErrVar <= 0)
    Rcpp::stop("measErrVar must be a positive number.\n");
  if (!is_finite(x0) || x0 <= 0)
    Rcpp::stop("x0 must be a positive number.\n");
  if (!is_finite(nu) || nu <= 0)
    Rcpp::stop("nu must be a positive number.\n");

  uword nt = timePnt.n_elem;
  // find the mean function
  rowvec meanVec = Rcpp::as<rowvec>(meanFunc(Rcpp::wrap(timePnt)));
  chk_mat(meanVec, "meanVec", "double");
  if (meanVec.n_elem != nt)
    Rcpp::stop("The length of output of meanVec is not equal to the length of timePnt.");

  // find the variance function
  rowvec varVec = Rcpp::as<rowvec>(varFunc(Rcpp::wrap(timePnt)));
  chk_mat(varVec, "varVec", "double");
  if (any(varVec < 0))
    Rcpp::stop("The variance function at timePnt must be greater than 0.");
  if (varVec.n_elem != nt)
    Rcpp::stop("The length of output of varVec is not equal to the length of timePnt.");

  // generate the correlation sequence
  vec corr = Rcpp::as<vec>(corrGen(timePnt - min(timePnt), corrType, x0, nu));

  // get the correlation funtion
  mat corrMat = zeros<mat>(nt, nt);
  for (uword i = 0; i < nt; i++)
    corrMat.col(i).tail(nt - i) += corr.head(nt - i);
  corrMat = symmatl(corrMat);

  // get the eigen values and eigen functions
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, corrMat);
  // get number of functional principal components
  uvec numFPC = find(cumsum(sort(eigval, "descend")) / sum(eigval) > 0.9, 1, "first");

  // generate random functional data
  mat loeveSum = randn(n, numFPC(0)) * eigvec.tail_cols(numFPC(0)).t();
  loeveSum.each_row() %= sqrt(varVec);
  mat errorTerm = sqrt(measErrVar) * randn(n, nt);
  mat funcData = loeveSum + errorTerm;
  funcData.each_row() += meanVec;

  return Rcpp::DataFrame::create(Rcpp::_["sampleID"] = repmat(linspace<uvec>(1, n, n), nt, 1),
                                 Rcpp::_["t"] = vectorise(repmat(timePnt, 1, n), 1).t(),
                                 Rcpp::_["y"] = vectorise(funcData));
}
