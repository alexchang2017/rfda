context("common functions")

context("1. factorial_f")
test_that("factorial_f", {
  expect_equal(rfda:::factorial_f(0L), factorial(0))
  expect_equal(rfda:::factorial_f(1L), factorial(1))
  expect_equal(rfda:::factorial_f(2L), factorial(2))
  expect_equal(rfda:::factorial_f(3L), factorial(3))
  expect_equal(rfda:::factorial_f(5L), factorial(5))
  expect_true(is.nan(rfda:::factorial_f(-1)))
  expect_error(rfda:::factorial_f(0.5))
  expect_error(rfda:::factorial_f(2.5))
})

context("2. unique_rows")
M <- matrix(c(1,1,2,2,3,3,1,1,1,1,4,5), 6, 2)
resM <- matrix(c(1,2,3,3,1,1,4,5), 4, 2)
M2 <- matrix(sample(1:10, 120, TRUE), 20, 6)
resM2 <- unique(M2)
resM2 <- resM2[do.call(order, lapply(1:ncol(resM2), function(i) resM2[, i])), ]
test_that("unique_rows", {
  expect_equal(unique_rows(M), resM)
  expect_equal(unique_rows(resM), resM)
  expect_equal(unique_rows(M2), resM2)
  expect_equal(unique_rows(matrix(M[1,], 1)), matrix(M[1,], 1))
  expect_error(unique_rows(rbind(M, c(1, NA))))
  expect_error(unique_rows(rbind(M, c(1, NaN))))
  expect_error(unique_rows(rbind(M, c(1, Inf))))
})

context("3. trapz")
xi <- 1:5
x <- c(1, 4, 9, 16, 25)
x2 <- matrix(c(1,4,9,16,25,1,8,27,64,125), 5)
test_that("trapz", {
  expect_equal(as.vector(trapz(x)), 42)
  expect_equal(as.vector(trapz(xi, x)), 42)
  expect_equal(as.vector(trapz(x2)), c(42, 162))
  expect_equal(as.vector(trapz(xi, x2)), c(42, 162))
  expect_error(trapz(c(NA, x)))
  expect_error(trapz(c(NaN, x)))
  expect_error(trapz(c(Inf, x)))
  expect_error(trapz(c(0, xi), x))
  expect_error(trapz(xi, c(0, x)))
  expect_error(trapz(c(1)))
  expect_error(trapz(c(0, xi), c(NA, x)))
  expect_error(trapz(c(0, xi), c(NaN, x)))
  expect_error(trapz(c(0, xi), c(Inf, x)))
  expect_error(trapz(c(NA, xi), c(0, x)))
  expect_error(trapz(c(NaN, xi), c(0, x)))
  expect_error(trapz(c(Inf, xi), c(0, x)))
})

context("4. binning data")
data2bin <- data.table(subId = rep(1:5, each = 5), timePnt = rep(1:5, 5),
                       variable = "y", value = 1:25)
binnedValue <- do.call(c, lapply(1:5, function(x) (x-1)*5 + seq(1.5, 4.5, 1.5)))
resDataDT <- data.table(subId = rep(1:5, each = 3), variable = "y",
                        value = binnedValue, timePnt = rep(c(5/3, 3, 13/3), 5))
test_that("binning data", {
  expect_equal(rfda:::binData(data2bin, 3), resDataDT)
  expect_error(rfda:::binData(data2bin, 0))
  expect_error(rfda:::binData(data2bin, -1))
  expect_error(rfda:::binData(data2bin, NA))
  expect_error(rfda:::binData(data2bin, NaN))
  expect_error(rfda:::binData(data2bin, Inf))
})

context("5. quantile")
x <- c(1, 4, 3, 7, 9 ,11, 15, 8, 2, 4, 8, 6)
x2 <- rnorm(100)
x3 <- sample(1:100, 100, TRUE)
p <- seq(0.05, 0.95, by = 0.05)
test_that("binning data", {
  expect_equal(quantile(x, p, type = 7), as.vector(rfda:::quantileCpp(x, p)), check.attributes = FALSE)
  expect_equal(quantile(x2, p, type = 7), as.vector(rfda:::quantileCpp(x2, p)), check.attributes = FALSE)
  expect_equal(quantile(x3, p, type = 7), as.vector(rfda:::quantileCpp(x3, p)), check.attributes = FALSE)
})

