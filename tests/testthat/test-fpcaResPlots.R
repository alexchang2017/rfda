context("test - fpcaResPlots")

data("irregularExData", package = 'rfda')
set.seed(100)
randCov <- matrix(rnorm(100), 10, 10)
randCov <- (randCov + t(randCov)) / 2
diag(randCov) <- diag(randCov) + 2
eigRes <- eigen(randCov, TRUE)
idx <- eigRes$values > 0
testFpcaRes <- list(data = irregularExData, id.var = "sampleID", time.var = "t",
                    eigVals = eigRes$values, numFPC = 4)
class(testFpcaRes) <- "fpcaRes"
test_that("test - fpcaResPlots", {
  expect_is(designPlot(irregularExData, "sampleID", "t"), "trellis")
  expect_is(designPlot(testFpcaRes), "trellis")
  expect_is(screePlot(eigRes$values[idx], sum(idx)), "trellis")
  expect_is(screePlot(eigRes$values[idx], 4, ylim = c(0, 0.8),
                      col.line = "red", pch.line = 17, lty.line = 4), "trellis")
  expect_is(screePlot(testFpcaRes), "trellis")
})
