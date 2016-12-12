context("test - FPCA")

data("sparseExData", package = 'rfda')
data("irregularExData", package = 'rfda')
data("regularExData", package = 'rfda')

context("1. test - checkSparsity")
test_that("checkSparsity", {
  expect_equal(checkSparsity(sparseExData, "sampleID", "t"), 0)
  expect_equal(checkSparsity(irregularExData, "sampleID", "t"), 1)
  expect_equal(checkSparsity(regularExData, "sampleID", "t"), 2)
})

context("2. test - FPCA check results (uni)")
sparsity <- checkSparsity(regularExData, "sampleID", "t")
bwCand <- bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1)
bwOpt <- gcvLocPoly1d(bwCand, regularExData$t, regularExData$y, rep(1, nrow(regularExData)), "gauss", 0, 1)
meanFunc <- locPoly1d(bwOpt, regularExData$t, regularExData$y, rep(1, nrow(regularExData)),
          sort(unique(regularExData$t)), "gauss", 0, 1)
udMF <- data.table(timePnt = sort(unique(regularExData$t)), value = meanFunc, variable = "y")
regularExData_multiVar <- regularExData %>>% data.table %>>% `[`(j = `:=`(y2 = y*sin(t)*3, y3 = y*cos(t)*2))
udCF <- list(
  "y-y" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y-y2" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y-y3" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y2-y2" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y2-y3" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y3-y3" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10))
)
sparseExDataTestError1 <- subset(sparseExData, sparseExData$sampleID == 1)
sparseExDataTestError2 <- subset(sparseExData, sparseExData$sampleID %in% c(1, 3))
test_that("test - FPCA check results (uni)", {
  expect_warning(FPCA(y ~ t, "sampleID", sparseExData), "Remove the observation with")
  expect_equal(suppressWarnings(FPCA(y ~ t, "sampleID", sparseExData)), 1)

  expect_equal(suppressWarnings(FPCA(y ~ t, "sampleID", sparseExData,
                                     list(weight = TRUE, ncpus = 1, numBins = 10))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", irregularExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData,
                    list(bwMean = data.table(variable = "y", value = 0.676065))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(bwMean = data.table(variable = "y", value = -2),
                                                           bwCov = data.table(variable1 = "y", variable2 = "y",
                                                                              value1 = -1, value2 = -1))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(userMeanFunc = udMF, userCovFunc = udCF[1],
                                                           bwCov = data.table(variable1 = "y", variable2 = "y",
                                                                              value1 = 0.5, value2 = 0.5))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(weight = TRUE, ncpus = 1, numBins = 10)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numBins = -1)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodNorm = "smoothCov")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodNorm = "quantile")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(newdata = 1:9)), 1)

  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(outPercent = 0.1)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(outPercent = 0.1, newdata = 1:9)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numFPC = 5)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numFPC = "AIC")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numFPC = "BIC")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodFPCS = "CE")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodFPCS = "IN")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodFPCS = "LS")), 1)
  expect_warning(FPCA(y ~ t, "sampleID", regularExData, list(numFPC = "AIC_R")), "so reset numFPC to")
})

context("2. test - FPCA check results (multi)")
test_that("test - FPCA check results (multi)", {
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(methodNorm = "smoothCov")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(methodNorm = "no")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(newdata = 1:9)), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(outPercent = 0.1)), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(outPercent = 0.1, newdata = 1:9)), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(numFPC = 5)), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(numFPC = "AIC")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(numFPC = "BIC")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(methodFPCS = "CE")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(methodFPCS = "IN")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(methodFPCS = "LS")), 1)
  expect_warning(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(numFPC = "AIC_R")), "so reset numFPC to")
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar,
                    list(bwMean = data.table(variable = c("y2", "y"), value = c(0.676065, -2)),
                         bwCov = data.table(variable1 = c("y", "y2", "y", "y", "y2"),
                                            variable2 = c("y", "y2", "y2", "y3", "y3"),
                                            value1 = c(-1, 0.5, -1, 0.7, -3),
                                            value2 = c(-1, 0.5, -1, 0.9, -3)))), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar,
                    list(userCovFunc = udCF,
                         bwCov = data.table(variable1 = c("y", "y2", "y3"), variable2 = c("y", "y2", "y3"),
                                            value1 = c(0.5, 0.5, 0.5), value2 = c(0.5, 0.5, 0.5)))), 1)
})

context("3. test - FPCA validate input")
test_that("test - FPCA validate input", {
  expect_error(FPCA(y ~ t, "sampleID", sparseExDataTestError1), "The number of curves cannot be less than 2!")
  expect_error(suppressWarnings(FPCA(y ~ t, "sampleID", sparseExDataTestError2)),
               "The number of curves cannot be less than 2!")
  expect_error(suppressWarnings(FPCA(y ~ t, "sampleID", sparseExData, list(methodFPCS = "IN"))),
               "The case methodFPCS = 'IN'")
  expect_error(FPCA(y ~ t, "sampleID", regularExData, list(userMeanFunc = udMF, userCovFunc = udCF[1])))
  expect_error(FPCA(y ~ t, "sampleID", regularExData, list(newdata = -1:2)), "The value of newdata")
  expect_error(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar,
                    list(userCovFunc = udCF,
                         bwCov = data.table(variable1 = c("y", "y2"), variable2 = c("y", "y2"),
                                            value1 = c(0.5, 0.5), value2 = c(0.5, 0.5)))),
               "The bandwidths of covariance function")
})
