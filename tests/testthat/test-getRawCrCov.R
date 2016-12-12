context("test - getRawCrCov")

context("1. test - check results")
dataDT <- data.table(variable = rep(1, 6), timePnt = rep(1:2, 3), subId = rep(1:3, each = 2), value.demean = -3:2)
test_that("test - validate input", {
  expect_equal(getRawCrCov(dataDT), data.table(variable1 = rep(1, 4), variable2 = rep(1, 4), t1 = rep(1:2, each = 2),
                                               t2 = rep(1:2, 2), sse = c(11, 8, 8, 8), cnt = rep(3, 4), weight = rep(1, 4)))
})

context("2. test - validate input")
dataDT1 <- data.table::copy(dataDT) %>>% data.table::setnames("variable", "variable1")
dataDT2 <- data.table::copy(dataDT)
data.table::set(dataDT2, 1L, 1L, NA)
test_that("test - validate input", {
  expect_error(getRawCrCov(dataDT1))
  expect_error(getRawCrCov(dataDT2))
})
