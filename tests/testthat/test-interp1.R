context("1D interpolation")

context("1. test interp1 - case 1")
x <- c(0.8, 0.3, 0.1, 0.6)
y <- x ** 2
ys <- y[order(x)]
xs <- x[order(x)]
xi <- seq(0.1, 0.8, len=5)
test_that("interp1 - find the 1D interpolation - case 1", {
  expect_equal(interp1(xs, ys, xi, 'linear'),
               c(0.01, 0.08, 0.225, 0.395, 0.64), tolerance = 1e-6)
  expect_equal(interp1(xs, ys, xi, 'spline'),
               c(0.01, 0.075625, 0.2025, 0.390625, 0.64), tolerance = 1e-6)
  expect_equal(interp1(c(0.1, 0.9), c(0.3, 0.8), c(0.1, 0.5, 0.9), 'spline'),
               c(0.3, 0.55, 0.8), tolerance = 1e-6)
})

context("2. test interp1 - unsorted data")
test_that("interp1 - test unsorted condition", {
  expect_equal(interp1(x, y, xi, 'linear'), interp1(xs, ys, xi, 'linear'))
  expect_equal(interp1(x, y, xi, 'spline'), interp1(xs, ys, xi, 'spline'))
  expect_message(interp1(x, y, xi, 'linear'), "x will be sorted.")
})

context("3. test interp1 - duplicated data")
test_that("interp1 - test duplicated data", {
  expect_equal(interp1(c(0.1, xs), c(0.05, ys), xi, 'spline'),
               interp1(xs, c(0.03, ys[2:4]), xi, 'spline'))
  expect_message(interp1(c(0.1, xs), as.matrix(c(0.05, ys)), xi, 'linear'),
                 "The values of y for duplicated values of x will be averaged.")
})


context("4. test spline_f")
test_that("interp1 - spline_f - case 1", {
  expect_equal(spline_f(xs, t(ys), xi),
               spline_f(xs, ys, xi), tolerance = 1e-6)
  expect_equal(interp1(xs, ys, xi, 'spline'),
               spline_f(xs, ys, xi), tolerance = 1e-6)
  expect_equal(spline_f(c(0.1, 0.9), c(0.3, 0.5, 0.7, 0.8), c(0.1, 0.5, 0.9)),
               c(0.5, 0.55, 0.7), tolerance = 1e-6)
  expect_equal(spline_f(c(0.1, 0.9), t(c(0.3, 0.5, 0.7, 0.8)), c(0.1, 0.5, 0.9)),
               c(0.5, 0.55, 0.7), tolerance = 1e-6)
  expect_equal(spline_f(xs, c(-0.02, ys, 0.81), xi),
               c(0.01, 0.07465259, 0.19448438, 0.39724731, 0.64), tolerance = 1e-6)
  expect_message(spline_f(c(0.1, xs), c(0.05, ys), xi),
                "The values of y for duplicated values of x will be averaged.")
  expect_message(as.vector(spline_f(x, as.matrix(y), xi)), "x will be sorted.")
})

context("5. test interp1 - input validation")
test_that("interp1 - test input error - case 1", {
  expect_error(interp1(c(0, xs), c(NA, ys), xi, 'spline'))
  expect_error(interp1(c(0, xs), c(NaN, ys), xi, 'spline'))
  expect_error(interp1(c(0, xs), c(Inf, ys), xi, 'spline'))
  expect_error(interp1(c(NA, xs), c(1, ys), xi, 'spline'))
  expect_error(interp1(c(NaN, xs), c(1, ys), xi, 'spline'))
  expect_error(interp1(c(Inf, xs), c(1, ys), xi, 'spline'))
  expect_error(interp1(xs, ys, c(NA, xi), 'spline'))
  expect_error(interp1(xs, ys, c(NaN, xi), 'spline'))
  expect_error(interp1(xs, ys, c(Inf, xi), 'spline'))
  expect_error(interp1(xs, c(ys, 1), xi, 'spline'))
  expect_error(interp1(c(0, xs), ys, xi, 'spline'))
  expect_error(interp1(xs[1], ys[1], xs[1], 'linear'))
  expect_error(interp1(xs[1], ys[1], xs[1], 'spline'))
  expect_error(interp1(xs, ys, xi, 'axspline'))
  expect_error(spline_f(xs, ys, c(NA, xi)))
  expect_error(spline_f(xs, ys, c(NaN, xi)))
  expect_error(spline_f(xs, ys, c(Inf, xi)))
  expect_error(spline_f(c(0, xs), c(NA, ys), xi))
  expect_error(spline_f(c(0, xs), c(NaN, ys), xi))
  expect_error(spline_f(c(0, xs), c(Inf, ys), xi))
  expect_error(spline_f(c(NA, xs), c(1, ys), xi))
  expect_error(spline_f(c(NaN, xs), c(1, ys), xi))
  expect_error(spline_f(c(Inf, xs), c(1, ys), xi))
  expect_error(spline_f(xs, c(ys, 1), xi))
  expect_error(spline_f(c(0, xs), ys, xi))
})

context("6. test interp1 - case 2")
x <- c(0.8, 0.3, 0.1, 0.6, 0.9, 0.5, 0.2, 0.0, 0.7, 1.0, 0.4)
y <- x**2
ys <- y[order(x)]
xs <- x[order(x)]
xi <- seq(0, 1, len=21)
test_that("interp1 - find the 1D interpolation - case 2", {
  expect_equal(interp1(xs, ys, xi, 'linear'),
               c(0.0, 0.005, 0.01, 0.025, 0.04, 0.065, 0.09, 0.125, 0.16, 0.205,
                 0.25, 0.305, 0.36, 0.425, 0.49, 0.565, 0.64, 0.725, 0.81, 0.905, 1.0),
               tolerance = 1e-6)
  expect_equal(interp1(xs, ys, xi, 'spline'),
               c(0.0, 0.0025, 0.01, 0.0225, 0.04, 0.0625, 0.09, 0.1225, 0.16, 0.2025, 0.25,
                 0.3025, 0.36, 0.4225, 0.49, 0.5625, 0.64, 0.7225, 0.81, 0.9025, 1.0),
               tolerance = 1e-6)
})

context("7. test interp1 - case 3")
x <- seq(0, 2*pi, pi/4)
y <- sin(x)
xi <- seq(0, 2*pi, pi/16)
test_that("interp1 - find the 1D interpolation - case 3", {
  expect_equal(interp1(x, y, xi, 'linear'), approx(x, y, xi)$y)
})
