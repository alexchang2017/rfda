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
test_that("FPCA", {
  expect_equal(FPCA(y ~ t, "sampleID", sparseExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", irregularExData), 1)
  expect_equal(FPCA(y ~ t, "sampleID", regularExData), 1)
})
