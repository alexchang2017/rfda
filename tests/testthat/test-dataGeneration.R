context("test - data generation")

x <- seq(-10, 10, len = 51L)

context("1. test - corrGen")
test_that("test - corrGen", {
  expect_equal(length(corrGen(x, "BesselJ")), 51L)
  expect_true(all(!is.na(corrGen(x, "BesselJ")) & is.finite(corrGen(x, "BesselJ"))))
  expect_equal(length(corrGen(x, "Matern")), 51L)
  expect_true(all(!is.na(corrGen(x, "Matern")) & is.finite(corrGen(x, "Matern"))))
  expect_equal(length(corrGen(x, "rq")), 51L)
  expect_true(all(!is.na(corrGen(x, "rq")) & is.finite(corrGen(x, "rq"))))
})

test_that("test - corrGen validate input", {
  expect_error(corrGen(c(x, NA), "BesselJ"))
  expect_error(corrGen(c(x, NaN), "BesselJ"))
  expect_error(corrGen(c(x, Inf), "BesselJ"))
  expect_error(corrGen(x, "Besselj"))
  expect_error(corrGen(x, "BesselJ", 0))
  expect_error(corrGen(x, "BesselJ", -1))
  expect_error(corrGen(x, "BesselJ", NA))
  expect_error(corrGen(x, "BesselJ", NaN))
  expect_error(corrGen(x, "BesselJ", Inf))
  expect_error(corrGen(x, "Matern", nu = 0))
  expect_error(corrGen(x, "Matern", nu = -1))
  expect_error(corrGen(x, "Matern", nu = NA))
  expect_error(corrGen(x, "Matern", nu = NaN))
  expect_error(corrGen(x, "Matern", nu = Inf))
})

context("2. test - funcDataGen")
n <- 10
test_that("test - funcDataGen", {
  expect_true(nrow(funcDataGen(n, x, function(x) sin(x),
                               function(x) rep(1, length(x)), "BesselJ")) == n * length(x))
  expect_true(nrow(funcDataGen(n, x, function(x) sin(x),
                               function(x) rep(1, length(x)), "Matern")) == n * length(x))
  expect_true(nrow(funcDataGen(n, x, function(x) sin(x),
                               function(x) rep(1, length(x)), "rq")) == n * length(x))
})

test_that("test - funcDataGen validate input", {
  expect_error(funcDataGen(n, c(x, NA), function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ"))
  expect_error(funcDataGen(n, c(x, NaN), function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ"))
  expect_error(funcDataGen(n, c(x, Inf), function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ"))
  expect_error(funcDataGen(n, x, function(x) 1,
                           function(x) rep(1, length(x)), "BesselJ"))
  expect_error(funcDataGen(n, x, function(x) 1/x,
                           function(x) rep(1, length(x)), "BesselJ"))
  expect_error(funcDataGen(n, x, function(x) rep(NA, length(x)),
                           function(x) rep(1, length(x)), "BesselJ"))
  expect_error(funcDataGen(n, x, function(x) sin(x), function(x) 1, "BesselJ"))
  expect_error(funcDataGen(n, x, function(x) sin(x), function(x) 1/x, "BesselJ"))
  expect_error(funcDataGen(n, x, function(x) sin(x), function(x) rep(NA, length(x)), "BesselJ"))
  expect_error(funcDataGen(n, x, function(x) sin(x), function(x) x - 5.5, "BesselJ"))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "Besselj"))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", -1))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", 0))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", NA))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", NaN))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", Inf))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", x0 = -1))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", x0 = 0))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", x0 = NA))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", x0 = NaN))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "BesselJ", x0 = Inf))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "Matern", nu = -1))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "Matern", nu = 0))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "Matern", nu = NA))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "Matern", nu = NaN))
  expect_error(funcDataGen(n, x, function(x) sin(x),
                           function(x) rep(1, length(x)), "Matern", nu = Inf))
})

context("3. test - sparsify")
DT <- funcDataGen(n, x, function(x) cos(x) + sin(x),
                  function(x) ((x - mean(x))**2)**(1/5)*2, "BesselJ")
test_that("test - funcDataGen", {
  expect_true(nrow(sparsify(DT, "sampleID", 0.1)) > 0)
  expect_true(nrow(sparsify(DT, "sampleID", 0.5)) > 0)
  expect_true(nrow(sparsify(DT, "sampleID", 0.9)) > 0)
  expect_true(nrow(sparsify(data.frame(aa = DT$sampleID), "aa", 0.2)) > 0)
})

test_that("test - funcDataGen validate input", {
  expect_error(sparsify(as.matrix(DT), "sampleID", 0.1))
  expect_error(sparsify(DT$sampleID, "sampleID", 0.1))
  expect_error(sparsify(DT, "sampleid", 0.1))
  expect_error(sparsify(DT, "sampleID", 0))
  expect_error(sparsify(DT, "sampleID", 1))
  expect_error(sparsify(DT, "sampleID", 1.5))
  expect_error(sparsify(DT, "sampleID", -1))
  expect_error(sparsify(DT, "sampleID", NA))
  expect_error(sparsify(DT, "sampleID", NaN))
  expect_error(sparsify(DT, "sampleID", Inf))
})
