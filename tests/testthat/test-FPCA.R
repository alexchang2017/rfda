context("FPCA")

data("sparseExData", package = 'rfda')
data("irregularExData", package = 'rfda')
data("regularExData", package = 'rfda')
regularExDataTestRemove <- data.table(regularExData)
regularExDataTestRemove <- rbind(regularExDataTestRemove[sampleID == 1][1, ],
                                 regularExDataTestRemove[sampleID == 2][1, ],
                                 regularExDataTestRemove[!sampleID %in% c(1, 2)])

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
regularExData_multiVar <- regularExData %>>% data.table %>>% `[`(j = `:=`(y2 = y*sin(t)*3, y3 = y*cos(t)*2))
udCF <- list(
  "y-y" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y-y2" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y-y3" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y2-y2" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y2-y3" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10)),
  "y3-y3" = matrix(rnorm(121), 11, 11, FALSE, list(0:10, 0:10))
)
test_that("FPCA", {
  expect_warning(FPCA(y ~ t, "sampleID", sparseExData), "Remove the observation with")
  expect_equal(suppressWarnings(FPCA(y ~ t, "sampleID", sparseExData)), 1)
  expect_equal(suppressWarnings(FPCA(y ~ t, "sampleID", sparseExData,
                                     list(weight = TRUE, ncpus = 1, numBins = 10))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", irregularExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData), 1)
  expect_warning(FPCA(y ~ t, "sampleID", regularExDataTestRemove), "Remove the observation with")
  expect_equal(FPCA(y ~ t, "sampleID", regularExData,
                    list(bwMean = data.table(variable = "y", value = 0.676065))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(bwMean = data.table(variable = "y", value = -2),
                                                           bwCov = data.table(variable1 = "y", variable2 = "y",
                                                                              value1 = -1, value2 = -1))), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(userMeanFunc = udMF, userCovFunc = udCF[1])), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(weight = TRUE, ncpus = 1, numBins = 10)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numBins = -1)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodNorm = "smoothCov")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(methodNorm = "quantile")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(newdata = 1:9)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(outPercent = 0.1)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(outPercent = 0.1, newdata = 1:9)), 1)
  expect_warning(FPCA(y ~ t, "sampleID", regularExData, list(outPercent = 0.3)), "Reset 'outPercent' to 25 percent!")
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numFPC = 5)), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numFPC = "AIC")), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData, list(numFPC = "BIC")), 1)
  expect_warning(FPCA(y ~ t, "sampleID", regularExData, list(numFPC = "AIC_R")), "so reset numFPC to")
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(methodNorm = "smoothCov")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(methodNorm = "no")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(newdata = 1:9)), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(outPercent = 0.1)), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(outPercent = 0.1, newdata = 1:9)), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(numFPC = 5)), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(numFPC = "AIC")), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(numFPC = "BIC")), 1)
  expect_warning(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar, list(numFPC = "AIC_R")), "so reset numFPC to")
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar,
                    list(bwMean = data.table(variable = c("y2", "y"), value = c(0.676065, -2)),
                         bwCov = data.table(variable1 = c("y", "y2", "y", "y", "y2"),
                                            variable2 = c("y", "y2", "y2", "y3", "y3"),
                                            value1 = c(-1, 0.5, -1, 0.7, -3),
                                            value2 = c(-1, 0.5, -1, 0.9, -3)))), 1)
  expect_equal(FPCA(y + y2 + y3 ~ t, "sampleID", regularExData_multiVar,
                    list(userCovFunc = udCF)), 1)
})
