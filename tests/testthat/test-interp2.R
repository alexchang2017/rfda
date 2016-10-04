context("2D interpolation")

load('interp2_res.rda') # results from MatLab

context("1. test interp2 - case 1")
A <- matrix(c(13, 5, 1, -1, 4, 6, 12, 3, 2), 3)
x <- c(0, 1, 4)
y <- c(10, 11, 12)
xi <- seq(min(x), max(x), len = 17)
yi <- seq(min(y), max(y), len = 26)
x_us <- c(0, 4, 1)
y_us <- c(10, 12, 11)
x2 <- x[c(1, 1:3)]
A2 <- A[ , c(1, 1:3)]
x3 <- x[c(1:3, 1)]
y2 <- y[c(1, 1:3)]
A4 <- A[c(1, 1:3), ]
test_that("interp2 - find the 2D interpolation - case 1", {
  expect_equal(interp2(x, y, A, xi, yi, 'linear'), V_1_l)
  expect_equal(interp2(x, y, A, xi, yi, 'spline'), V_1_s)
  expect_equal(interp2(x2, y, A2, xi, yi, 'linear'), V_1_l)
  expect_equal(interp2(x2, y, A2, xi, yi, 'spline'), V_1_s)

})

context("2. test interp2 - unsorted data")
test_that("interp2 - test unsorted condition", {
  expect_equal(interp2(x_us, y_us, A[order(y_us), order(x_us)], xi, yi, 'linear'), V_1_l)
  expect_equal(interp2(x_us, y_us, A[order(y_us), order(x_us)], xi, yi, 'spline'), V_1_s)
  expect_message(interp2(x_us, y, A[ , order(x_us)], xi, yi, 'linear'), "x will be sorted.")
  expect_message(interp2(x, y_us, A[ , order(x_us)], xi, yi, 'linear'), "y will be sorted.")
})

context("3. test interp2 - duplicated data")
test_that("interp2 - test duplicated data", {
  expect_message(interp2(x2, y, A2, xi, yi, 'spline'),
                 "The values of v for duplicated values of x will be averaged.")
  expect_message(interp2(x, y2, A4, xi, yi, 'spline'),
                 "The values of v for duplicated values of y will be averaged.")
  expect_equal(interp2(x, y2, A4, xi, yi, 'spline'), V_1_s)
  expect_equal(interp2(x3, y, A2[, rank(x3, ties.method="first")], xi, yi, 'linear'), V_1_l)
})

context("4. test interp2 - case 2")
x <- seq(-3, 3, 1)
xm <- expand.grid(x, x)
z <- 3*(1-xm[,1])^2.*exp(-(xm[,1]^2) - (xm[,2]+1)^2) -
  10*(xm[,1]/5 - xm[,1]^3 - xm[,2]^5)*exp(-xm[,1]^2-xm[,2]^2) -
  1/3*exp(-(xm[,1]+1)^2 - xm[,2]^2)
Z <- matrix(z, length(x), byrow = TRUE)
xi <- seq(-3, 3, 0.25)
test_that("interp2 - find the 2D interpolation - case 2", {
  expect_equal(interp2(x, x, Z, xi, xi, 'linear'), V_2_l)
  expect_equal(interp2(x, x, Z, xi, xi, 'spline'), V_2_s)
})

context("6. test interp2 - case 3")
x <- c(0.1, 0.2, 0.8)
y <- c(0.1, 0.3, 0.6, 0.8)
v <- matrix(c(0.2,0.3,0.4,0.5,0.5,0.4,0.4,0.4,0.2, 0.4, 0.5, 0.4), 4)
xi <- seq(0.1, 0.8, len = 8)
test_that("interp2 - find the 2D interpolation - case 3", {
  expect_equal(interp2(x, y, v, xi, xi, 'linear'), V_3_l)
  expect_equal(interp2(x, y, v, xi, xi, 'spline'), V_3_s)
})

context("7. test interp2 - input validation")
v21 <- v
v21[1] <- NA
v22 <- v
v22[1] <- NaN
v23 <- v
v23[1] <- Inf
test_that("interp2 - test input error - case 3", {
  expect_error(interp2(c(0, x), y, v, xi, xi, 'spline'))
  expect_error(interp2(c(NA, x[1:2]), y, v, xi, xi, 'spline'))
  expect_error(interp2(c(NaN, x[1:2]), y, v, xi, xi, 'spline'))
  expect_error(interp2(c(Inf, x[1:2]), y, v, xi, xi, 'spline'))
  expect_error(interp2(x, c(0, y), v, xi, xi, 'spline'))
  expect_error(interp2(x, c(NA, y[1:2]), v, xi, xi, 'spline'))
  expect_error(interp2(x, c(NaN, y[1:2]), v, xi, xi, 'spline'))
  expect_error(interp2(x, c(Inf, y[1:2]), v, xi, xi, 'spline'))
  expect_error(interp2(x, y, v21, xi, xi, 'spline'))
  expect_error(interp2(x, y, v22, xi, xi, 'spline'))
  expect_error(interp2(x, y, v23, xi, xi, 'spline'))
  expect_error(interp2(x, y, v, c(NA, xi), xi, 'spline'))
  expect_error(interp2(x, y, v, c(NaN, xi), xi, 'spline'))
  expect_error(interp2(x, y, v, c(Inf, xi), xi, 'spline'))
  expect_error(interp2(x, y, v, xi, c(NA, xi), 'spline'))
  expect_error(interp2(x, y, v, xi, c(NaN, xi), 'spline'))
  expect_error(interp2(x, y, v, xi, c(Inf, xi), 'spline'))
  expect_error(interp2(x, y, v, xi, xi, 'splinea'))
})

