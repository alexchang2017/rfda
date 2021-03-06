% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{locLinear2d}
\alias{locLinear2d}
\title{Two-dimensional kernel linear local smoother}
\usage{
locLinear2d(bandwidth, x, y, w, count, out1, out2, kernel)
}
\arguments{
\item{bandwidth}{A numeric vector with two values. The kernel smoothing parameters.}

\item{x}{A matrix, the variable of of x-axis and y-axis.}

\item{y}{A vector, the variable of of z-axis. \code{y[i]} is corresponding value of \code{x[i, ]}.}

\item{w}{A vector, the weight of data. \code{w[i]} is corresponding value of \code{x[i,]}.}

\item{count}{A vector, the number of observations at \code{x[i, ]}.}

\item{out1}{A vector, the output grid of x-coordinate. It should be a sorted vecotr.}

\item{out2}{A vector, the output grid of y-coordinate. It should be a sorted vecotr.}

\item{kernel}{A string. It could be 'gauss', 'gaussvar', 'epan' or 'quar'.}
}
\value{
A smoothed covariance estimated by two-dimensional kernel local linear smoother.
}
\description{
Perform two-dimensional kernel linear local smoother for data \code{(x,y)} with weight \code{w} on \code{xout}.
}
\examples{
data("regularExData", package = 'rfda')
sparsity <- checkSparsity(regularExData, "sampleID", "t")
bwCand <- bwCandChooser(regularExData, "sampleID", "t", sparsity, "gauss", 1)
w <- rep(1, nrow(regularExData))
bwOpt <- gcvLocPoly1d(bwCand, regularExData$t, regularExData$y, w, "gauss", 0, 1)
bwOpt <- adjGcvBw(bwOpt, sparsity, "gauss", 0)
xout <- sort(unique(regularExData$t))
meanFunc <- locPoly1d(bwOpt, regularExData$t, regularExData$y, w, xout, "gauss", 0, 1)
require(data.table)
require(pipeR)
setDT(regularExData, key = "t")
dataDT <- merge(regularExData, data.table(mf = meanFunc, t = xout, key = "t"), by = "t") \%>>\%
  `[`( , `:=`(y = y - mf, variable = "y")) \%>>\%
  setnames(c("t", "y", "sampleID"), c("timePnt", "value.demean", "subId"))
RawCov <- getRawCrCov(dataDT)

xout2 <- seq(min(regularExData$t), max(regularExData$t), len = 30)
RawCovNoDiag <- RawCov[t1 != t2]
covFunc <- locLinear2d(c(1, 1), as.matrix(RawCovNoDiag[ , .(t1, t2)]), RawCovNoDiag$sse,
  RawCovNoDiag$weight, RawCovNoDiag$cnt, xout2, xout2, "gauss")
}

