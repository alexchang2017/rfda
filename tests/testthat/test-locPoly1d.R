context("test - locPoly1d")

context("test - bwCandChooser")
data("regularExData", package = 'rfda')
data("irregularExData", package = 'rfda')
data("sparseExData", package = 'rfda')
test_that("lwls - kernel", {
  expect_equal(bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1),
               c(0.789474, 0.897346, 1.019959, 1.159325, 1.317734, 1.497787, 1.702443,
                 1.935063, 2.199467, 2.500000), tolerance = 1e-6)
  expect_equal(bwCandChooser(irregularExData, "sampleID", "t", 1, "gauss", 1),
               c(1.052632, 1.158822, 1.275725, 1.404422, 1.546101, 1.702074, 1.873781,
                 2.062810, 2.270908, 2.500000), tolerance = 1e-6)
  expect_equal(bwCandChooser(sparseExData, "sampleID", "t", 0, "gauss", 1),
               c(0.264199, 0.338967, 0.434895, 0.557970, 0.715876, 0.918469, 1.178395,
                 1.511881, 1.939743, 2.488689), tolerance = 1e-6)
})

context("test - locPoly1d with same weights")
bw <- 0.25
x <- seq(0, 1, 0.1)
y <- x ** 2 *2 + 3 * x
w <- rep(1, 11)
test_that("locPoly1d - with same weights", {
  expect_equal(as.vector(locPoly1d(bw, x, y, w, x, "epan", 0, 1)),
               c(-0.004684, 0.337088, 0.706824, 1.106824, 1.546824, 2.026824, 2.546824,
                 3.106824, 3.706824, 4.337088, 4.995316), tolerance = 1e-6)
  expect_equal(as.vector(locPoly1d(bw, x, y, w, x, "quar", 0, 1)),
               c(-0.002781, 0.334288, 0.698334, 1.098334, 1.538334, 2.018334, 2.538334,
                 3.098334, 3.698334, 4.334288, 4.997219), tolerance = 1e-6)
  expect_equal(as.vector(locPoly1d(bw, x, y, w, x, "gauss", 0, 1)),
               c(-0.059163, 0.325037, 0.731015, 1.161760, 1.619430, 2.105204, 2.619430,
                 3.161760, 3.731015, 4.325037, 4.940837), tolerance = 1e-6)
  expect_equal(as.vector(locPoly1d(bw, x, y, w, x, "gaussvar", 0, 1)),
               c(0.006464, 0.345469, 0.725682, 1.140981, 1.590669, 2.074075, 2.590669,
                 3.140981, 3.725682, 4.345469, 5.006464), tolerance = 1e-6)
})

context("test - locPoly1d validate inputs")
test_that("test - locPoly1d validate inputs", {
  expect_true(is.na(locPoly1d(0.01, x, y, w, x, "epan", 0, 1)))
  expect_error(locPoly1d(bw, x, y, w, x, "guass", 0, 1), 'Unsupported kernel')
  expect_error(locPoly1d(-0.2, x, y, w, x, "gauss", 0, 1))
  expect_error(locPoly1d(NA, x, y, w, x, "gauss", 0, 1))
  expect_error(locPoly1d(NaN, x, y, w, x, "gauss", 0, 1))
  expect_error(locPoly1d(Inf, x, y, w, x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(1, x), y, w, x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(NA, x), c(1, y), c(1, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(NaN, x), c(1, y), c(1, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(Inf, x), c(1, y), c(1, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, x, c(1, y), w, x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(1, x), c(NA, y), c(1, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(1, x), c(NaN, y), c(1, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(1, x), c(Inf, y), c(1, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, x, y, c(1, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(1, x), c(1, y), c(NA, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(1, x), c(1, y), c(NaN, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, c(1, x), c(1, y), c(Inf, w), x, "gauss", 0, 1))
  expect_error(locPoly1d(bw, x, y, w, c(NA, x), "gauss", 0, 1))
  expect_error(locPoly1d(bw, x, y, w, c(NaN, x), "gauss", 0, 1))
  expect_error(locPoly1d(bw, x, y, w, c(Inf, x), "gauss", 0, 1))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", -1, 1))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", NA, 1))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", NaN, 1))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", Inf, 1))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", 0, -1))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", 0, NA))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", 0, NaN))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", 0, Inf))
  expect_error(locPoly1d(bw, x, y, w, x, "gauss", 2, 1))
})

context("test - locPoly1d with different weights")
w <- c(rep(1, 9), 0, 0.5)
test_that("locPoly1d - with different weights", {
  expect_equal(as.vector(locPoly1d(bw, x, y, w, x, "epan", 0, 1)),
               c(-0.004684, 0.337088, 0.706824, 1.106824, 1.546824, 2.026824, 2.546824,
                 3.097088, 3.701979, 4.342951, 5.000000), tolerance = 1e-6)
  expect_equal(as.vector(locPoly1d(bw, x, y, w, x, "quar", 0, 1)),
               c(-0.002781, 0.334288, 0.698334, 1.098334, 1.538334, 2.018334, 2.538334,
                 3.094288, 3.691107, 4.341831, 5.000000), tolerance = 1e-6)
  expect_equal(as.vector(locPoly1d(bw, x, y, w, x, "gauss", 0, 1)),
               c(-0.057596, 0.326707, 0.731523, 1.158736, 1.609788, 2.086293, 2.590595,
                 3.125604, 3.693361, 4.292672, 4.917981), tolerance = 1e-6)
  expect_equal(as.vector(locPoly1d(bw, x, y, w, x, "gaussvar", 0, 1)),
               c(-0.002448, 0.341471, 0.725772, 1.142526, 1.588922, 2.064847, 2.573735,
                 3.122181, 3.716123, 4.355114, 5.038236), tolerance = 1e-6)
})

context("test - gcv_locPoly1d")

regBwCand <- bwCandChooser(regularExData, "sampleID", "t", 2, "gauss", 1)
irrBwCand <- bwCandChooser(irregularExData, "sampleID", "t", 1, "gauss", 1)
spsBwCand <- bwCandChooser(sparseExData, "sampleID", "t", 0, "gauss", 1)

regWeight <- rep(1, nrow(regularExData))
irrWeight <- rep(1, nrow(irregularExData))
spsWeight <- rep(1, nrow(sparseExData))

test_that("gcv_locPoly1d", {
  expect_equal(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, 1), 0.789474, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gaussvar", 0, 1), 0.789474, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "epan", 0, 1), 1.019959, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "quar", 0, 1), 1.019959, tolerance = 1e-6)

  expect_equal(gcv_locPoly1d(irrBwCand, irregularExData$t, irregularExData$y,
                             irrWeight, "gauss", 0, 1), 1.052632, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(irrBwCand, irregularExData$t, irregularExData$y,
                             irrWeight, "gaussvar", 0, 1), 1.052632, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(irrBwCand, irregularExData$t, irregularExData$y,
                             irrWeight, "epan", 0, 1), 1.052632, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(irrBwCand, irregularExData$t, irregularExData$y,
                             irrWeight, "quar", 0, 1), 1.052632, tolerance = 1e-6)

  expect_equal(gcv_locPoly1d(spsBwCand, sparseExData$t, sparseExData$y,
                             spsWeight, "gauss", 0, 1), 0.434895, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(spsBwCand, sparseExData$t, sparseExData$y,
                             spsWeight, "gaussvar", 0, 1), 0.557970, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(spsBwCand, sparseExData$t, sparseExData$y,
                             spsWeight, "epan", 0, 1), 0.918469, tolerance = 1e-6)
  expect_equal(gcv_locPoly1d(spsBwCand, sparseExData$t, sparseExData$y,
                             spsWeight, "quar", 0, 1), 1.178395, tolerance = 1e-6)
})

context("test - locPoly1d validate inputs")
test_that("test - locPoly1d validate inputs", {
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                               regWeight, "guass", 0, 1), 'Unsupported kernel')
  expect_error(gcv_locPoly1d(-0.2, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, 1))
  expect_error(gcv_locPoly1d(NA, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, 1))
  expect_error(gcv_locPoly1d(NaN, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, 1))
  expect_error(gcv_locPoly1d(Inf, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(1, regularExData$t), regularExData$y,
                             regWeight, "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(NA, regularExData$t), c(1, regularExData$y),
                             c(1, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(NaN, regularExData$t), c(1,regularExData$y),
                             c(1, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(Inf, regularExData$t), c(1,regularExData$y),
                             c(1, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, x, c(1, y), w, x, "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(1, regularExData$t), c(NA, regularExData$y),
                             c(1, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(1, regularExData$t), c(NaN,regularExData$y),
                             c(1, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(1, regularExData$t), c(Inf,regularExData$y),
                             c(1, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             c(1, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(1, regularExData$t), c(1, regularExData$y),
                             c(NA, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(1, regularExData$t), c(1, regularExData$y),
                             c(NaN, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, c(1, regularExData$t), c(1, regularExData$y),
                             c(Inf, regWeight), "gauss", 0, 1))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", -1, 1))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", NA, 1))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", NaN, 1))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", Inf, 1))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, -1))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, NA))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, NaN))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", 0, Inf))
  expect_error(gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "gauss", 2, 1))
})

context("test - adjGcvBw1d")

reg_bw_gauss <- gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                              regWeight, "gauss", 0, 1)
reg_bw_epan <- gcv_locPoly1d(regBwCand, regularExData$t, regularExData$y,
                             regWeight, "epan", 0, 1)
irr_bw_gauss <- gcv_locPoly1d(irrBwCand, irregularExData$t, irregularExData$y,
                              irrWeight, "gauss", 0, 1)
irr_bw_epan <- gcv_locPoly1d(irrBwCand, irregularExData$t, irregularExData$y,
                             irrWeight, "epan", 0, 1)
sps_bw_gauss <- gcv_locPoly1d(spsBwCand, sparseExData$t, sparseExData$y,
                              spsWeight, "gauss", 0, 1)
sps_bw_epan <- gcv_locPoly1d(spsBwCand, sparseExData$t, sparseExData$y,
                             spsWeight, "epan", 0, 1)

test_that("adjGcvBw1d", {
  expect_equal(adjGcvBw1d(reg_bw_gauss, 2, "gauss", 0), 0.868421, tolerance = 1e-6)
  expect_equal(adjGcvBw1d(reg_bw_epan, 2, "epan", 0), 1.121955, tolerance = 1e-6)
  expect_equal(adjGcvBw1d(irr_bw_gauss, 1, "gauss", 0), 1.157895, tolerance = 1e-6)
  expect_equal(adjGcvBw1d(irr_bw_epan, 1, "epan", 0), 1.157895, tolerance = 1e-6)
  expect_equal(adjGcvBw1d(sps_bw_gauss, 0, "gauss", 0), 0.478385, tolerance = 1e-6)
  expect_equal(adjGcvBw1d(sps_bw_epan, 0, "epan", 0), 1.010316, tolerance = 1e-6)
})

