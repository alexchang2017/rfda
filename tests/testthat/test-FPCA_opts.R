context("FPCA_opts")

context("1. get_FPCA_opts")
FPCA_opts_uni <- get_FPCA_opts(1)
FPCA_opts_mul <- get_FPCA_opts(2)
allOptNames <- c("bwMean", "bwCov", "bwNumGrid", "bwKernel", "numBins", "errTerm", "numGrid", "weight",
                 "numFPC", "FVE_threshold", "maxNumFPC", "methodFPCS", "shrink", "varErr", "outPercent",
                 "methodNorm", "quantile_probs", "ncpus", "userMeanFunc", "userCovFunc")
test_that("test get_FPCA_opts", {
  expect_equal(names(FPCA_opts_uni), allOptNames)
  expect_equal(names(FPCA_opts_mul), allOptNames)
  expect_equal(FPCA_opts_uni$methodNorm, "no")
  expect_equal(FPCA_opts_mul$methodNorm, "quantile")
})

context("2. chk_FPCA_opts")
default_FPCA_opts <- get_FPCA_opts(1)
testFPCAOpts <- function(x) modifyList(default_FPCA_opts, x) %>>% (rfda:::chk_FPCA_opts(.))
test_that("test get_FPCA_opts", {
  expect_warning(testFPCAOpts(list(numFPC = 54)), "Reset it to be numGrid-2 now!")
  expect_warning(testFPCAOpts(list(numFPC = "AIC", maxNumFPC  = 54)), "Reset it to be numGrid-2 now!")
  expect_warning(testFPCAOpts(list(shrink = TRUE, errTerm  = FALSE)), "Reset to shrink = FALSE now!")
  expect_error(testFPCAOpts(list(varErr = -0.5)))
})
