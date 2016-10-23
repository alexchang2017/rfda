context("FPCA")

data("sparseExData", package = 'rfda')
data("irregularExData", package = 'rfda')
data("regularExData", package = 'rfda')

context("1. checkSparsity")
test_that("checkSparsity", {
  expect_equal(checkSparsity(sparseExData, "sampleID", "t"), 0)
  expect_equal(checkSparsity(irregularExData, "sampleID", "t"), 1)
  expect_equal(checkSparsity(regularExData, "sampleID", "t"), 2)
})

context("2. FPCA")
sparsity <- checkSparsity(regularExData, "sampleID", "t")
bwCand <- bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1)
bwOpt <- gcvLocPoly1d(bwCand, regularExData$t, regularExData$y, rep(1, nrow(regularExData)), "gauss", 0, 1)
meanFunc <- locPoly1d(bwOpt, regularExData$t, regularExData$y, rep(1, nrow(regularExData)),
          sort(unique(regularExData$t)), "gauss", 0, 1)
udMF <- data.table(timePnt = sort(unique(regularExData$t)), value = meanFunc, variable = "y")
test_that("FPCA", {
  expect_equal(FPCA(y ~ t, "sampleID", sparseExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", sparseExData, list(weight = TRUE, npus = 1, numBins = 10)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", irregularExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(userMeanFuncList = udMF)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(weight = TRUE, npus = 1, numBins = 10)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodNorm = "smoothCov")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodNorm = "quantile")), 1)
})
