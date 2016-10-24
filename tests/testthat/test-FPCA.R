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
regularExData_multiVar <- regularExData %>>% data.table %>>% `[`(j = `:=`(y2 = y*sin(t), y3 = y*cos(t)))
test_that("FPCA", {
  expect_equal(FPCA(y ~ t, "sampleID", sparseExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", sparseExData, list(weight = TRUE, ncpus = 1, numBins = 10)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", irregularExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(bwMean = data.table(variable = "y", value = 0.676065))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(bwMean = data.table(variable = "y", value = -2))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(userMeanFunc = udMF)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(weight = TRUE, ncpus = 1, numBins = 10)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numBins = -1)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodNorm = "smoothCov")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodNorm = "quantile")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar,
                    list(bwMean = data.table(variable = c("y2", "y"), value = c(0.676065, -2)))), 1)
})
