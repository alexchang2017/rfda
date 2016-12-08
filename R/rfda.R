#' The rfda package: summary information
#'
#' to be done.
#'
#' @references
#' \enumerate{
#'   \item Yao, F., Muller, H.G., Clifford, A.J., Dueker, S.R., Follett, J., Lin, Y., Buchholz, B., Vogel, J.S. (2003). Shrinkage estimation for functional principal component scores, with application to the population kinetics of plasma folate. Biometrics 59, 676-685.
#'   \item Yao, F., Muller, H.G., Wang, J.L. (2005). Functional data analysis for sparse longitudinal data. J. American Statistical Association 100, 577-590.
#'   \item Yao, F., Muller, H.G., Wang, J.L. (2005). Functional Linear Regression Analysis for Longitudinal Data. Annals of Statistics 33, 2873-2903.
#'   \item Muller, H.G. (2005). Functional modeling and classification of longitudinal data. Scandinavian J. Statistics 32, 223-240.
#'   \item Muller, H.G., Stadtmuller, U., Yao, F. (2006). Functional variance processes. Journal of the American Statistical Association 101, 1007-1018.
#'   \item Leng, X., Muller, H.G. (2006). Classification using functional data analysis for temporal gene expression data. Bioinformatics 22, 68-76.
#'   \item Leng, X., Muller, H.G. (2006). Time ordering of gene co-expression. Biostatistics 7, 569-584.
#'   \item Chiou, J., Muller, H.G. (2007). Diagnostics for functional regression via residual processes. Computational Statistics & Data Analysis 51, 4849-4863.
#'   \item Muller, H.G. (2009). Functional modeling of longitudinal data. In: Longitudinal Data Analysis (Handbooks of Modern Statistical Methods), Ed. Fitzmaurice, G., Davidian, M., Verbeke, G., Molenberghs, G., Wiley, New York, 223--252.
#'   \item Hall, P., Muller, H.G., Yao, F. (2008). Modeling sparse generalized longitudinal observations via latent Gaussian processes. Journal of the Royal Statistical Society B 70, 703-723.
#'   \item Muller, H.G., Chiou, J.M., Leng, X. (2008). Inferring gene expression dynamics via functional regression analysis. BMC Bioinformatics 9:60.
#'   \item Peng, J., Muller, H.G. (2008). Distance-based clustering of sparsely observed stochastic processes, with applications to online auctions. Annals of Applied Statistics 2, 1056-1077.
#'   \item Tang, R., Muller, H.G. (2008). Pairwise curve synchronization for high-dimensional data. Biometrika 95, 875-889.
#'   \item Tang, R., Muller, H.G. (2009). Time-synchronized clustering of gene expression trajectories. Biostatistics 10, 32-45.
#'   \item Liu, B., Muller, H.G. (2009). Estimating derivatives for samples of sparsely observed functions, with application to on-line auction dynamics. J. American Statistical Association 104, 704-717.
#'   \item Muller, H.G., Yang, W. (2010). Dynamic relations for sparsely sampled Gaussian processes. Test 19, 1-29.
#'   \item Muller, H.G., Yao, F. (2010). Empirical dynamics for longitudinal data. Annals of Statistics 38, 3458?486.
#'   \item Chiou, J.M., Muller, H.G., Wang, J.L. (2003). Functional quasi-likelihood regression with smooth random effects. J. Royal Statistical Society B 65, 405-423.
#'   \item Chiou, J.M., Muller, H.G., Wang, J.L. (2004). Functional response models. Statistica Sinica 14, 675-693.
#'   \item Chiou, J., Muller, H.G. (2004). Quasi-likelihood regression with multiple indices and smooth link and variance functions. Scandinavian J. Statistics 31, 367-386.
#'   \item Yang, W., Muller, H.G., Stadtmuller, U. (2011). Functional singular component analysis. J. Royal Statistical Society B 73, 303?324.
#'   \item Dubin, J., Muller, H.G. (2005). Dynamical correlation for multivariate longitudinal data. J. American Statistical Association 100, 872-881.
#'   \item Chen, K., Muller, H.G. (2011). Conditional quantile analysis when covariates are functions, with application to growth data. J. Royal Statistical Society B 74, 67?9.
#'   \item Chen, K., Chen, K., Muller, H.G., Wang, J.L. (2011). Stringing high-dimensional data for functional analysis. J. American Statistical Association 106, 275-284.
#'   \item Chen, D., Muller, H.G. (2012). Nonlinear manifold representations for functional data. Annals of Statistics 40, 1-29.
#'   \item Muller, H.G., Sen, R., Stadtmuller, U. (2011). Functional Data Analysis for Volatility. J. Econometrics 165, 233-245.
#'   \item Jiang, C.R. and Wang, J.L. (2010): Covariate Adjusted Functional Principal Components Analysis for Longitudinal Data, The Annals of Statistics, 38, 1194-1226.
#'   \item Jiang, C.R. and Wang, J.L. (2011): Functional Single Index Model for Longitudinal Data, The Annals of Statistics, 39, 362-388.
#'   \item Chen, K. and Muller, H.G. (2012). Modeling Repeated Functional Observations, Journal of the American Statistical Association, 107, 1599-1609.
#'   \item Gottlieb, A. and Muller, H.G. (2012). A Stickiness Coefficient for Longitudinal Data, Computational Stat. and Data Analysis, i56, 4000-4010.
#'   \item Verzelen, N., Tao, W. and Muller, H.G. (2012). Inferring Stochastic Dynamics from Functional Data. Biometrika, 99, 533-550.
#'   \item Chiou, J.M. and Li, P.L. (2007). Functional clustering and identifying substructures of longitudinal data[J]. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 69, 679-699.
#'   \item Chiou, J.M., Chen, Y.T. and Yang, Y.F. (2014), Multivariate Functional Principal Component Analysis: A Normalization Approach. Statistica Sinica, 24, 1571-1596.
#'   \item Chiou, J.M., Yang, Y.F and Chen, Y.T.. (2016), Multivariate functional linear regression and prediction. Journal of Multivariate Analysis, 146, 301-312
#' }
#'
#' @docType package
#' @name rfda
#' @useDynLib rfda
#' @import assertthat
#' @importFrom Rcpp cppFunction sourceCpp
#' @importFrom pipeR %>>%
#' @importFrom data.table :=
#' @importFrom utils globalVariables data
#' @importFrom RcppParallel RcppParallelLibs
utils::globalVariables(c(".", "%>>%", ":=", "~<-"))
# variables used in data.table
utils::globalVariables(c("cnt", "idx_agg", "subId", "t1", "t2", "timePnt", "value", "value.mean",
                         "value.var1", "value.var2", "variable", "variable1", "variable2", "weight",
                         "value1", "value2", "value.demean", "value.var", "rn", "tmp"))

# Copyright (c) 2015, Hans-Georg Mueller and Jane-Ling Wang
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
# OF SUCH DAMAGE.
