context("test - locLinear2d")

load("covRes.rda") # results from MatLab
data("regularExData", package = "rfda")
data("irregularExData", package = "rfda")
data("sparseExData", package = "rfda")

testDataRegular <- data.table(y = c(seq(0.1, 0.9, 0.1), c(0.9,1,1)), t = rep(1:3, 4),
                              sampleID = c(rep(1:4, each=3)))
testDataIrregular <- data.table(y = seq(0.1, 1, 0.1), t = c(1,2,4,1:3,3,4,1,4),
                                sampleID = c(rep(1:2, each=3), rep(3:4, each=2)))
testDataSparse <- data.table(y = seq(0.1, 1, 0.1), t = c(1,2,4,3,5,7,6,10,8,9),
                             sampleID = c(rep(1:2, each=3), rep(3:4, each=2)))
testDT_list <- suppressMessages(lapply(list(testDataRegular, testDataIrregular, testDataSparse), function(dt){
  sparsity <- checkSparsity(dt, "sampleID", "t")
  bwCand <- bwCandChooser(dt, "sampleID", "t", sparsity, "gauss", 1)
  w <- rep(1, nrow(dt))
  bwOpt <- gcvLocPoly1d(bwCand, dt$t, dt$y, w, "gauss", 0, 1)
  bwOpt <- rfda:::adjGcvBw1d(bwOpt, sparsity, "gauss", 0)
  xout <- sort(unique(dt$t))
  meanFunc <- locPoly1d(bwOpt, dt$t, dt$y, w, xout, "gauss", 0, 1)

  merge(dt, data.table(mf = meanFunc, t = xout), by = "t") %>>%
    `[`( , `:=`(y = y - mf, variable = "y")) %>>%
    setnames(c("t", "y", "sampleID"), c("timePnt", "value", "subId")) %>>%
    `[`( , mf := NULL)
}))

compareCovResList <- lapply(list(rcov_case2, rcov_case1, rcov_case0), function(m){
  data.table(m) %>>% setnames(c("t1", "t2", "sse")) %>>% setorder(t1, t2)
})

context("1. test - getCrCov")
rawCovList <- lapply(testDT_list, rfda:::getRawCrCov)
test_that("test - bwCandChooser2", {
  expect_equal(rawCovList[[1]][ , .(t1,t2,sse)], compareCovResList[[1]] %>>% `[`( , sse := sse * 4))
  expect_equal(rawCovList[[2]][ , .(t1,t2,sse)], compareCovResList[[2]])
  expect_equal(rawCovList[[3]][ , .(t1,t2,sse)], compareCovResList[[3]])
})

context("2. test - bwCandChooser2")
allGridList <- lapply(list(regularExData, irregularExData, sparseExData), function(df){
  data.table(df) %>>% `[`( , .(t1 = rep(t, length(t)), t2 = rep(t, each=length(t))), by = .(sampleID))
})
test_that("test - bwCandChooser2", {
  expect_equal(bwCandChooser2(allGridList[[1]], "sampleID", c("t1", "t2"), 2, "gauss", 1)[ , 1],
               c(0.315789, 0.397407, 0.500119, 0.629378, 0.792045, 0.996754, 1.254371,
                 1.578570, 1.986561, 2.500000), tolerance = 1e-6)
  expect_equal(bwCandChooser2(allGridList[[2]], "sampleID", c("t1", "t2"), 1, "gauss", 1)[ , 1],
               c(0.210526, 0.277147, 0.364850, 0.480306, 0.632297, 0.832387, 1.095794,
                 1.442556, 1.899049, 2.500000), tolerance = 1e-6)
  expect_equal(bwCandChooser2(allGridList[[3]], "sampleID", c("t1", "t2"), 0, "gauss", 1)[ , 1],
               c(0.995476, 1.102164, 1.220286, 1.351068, 1.495867, 1.656184, 1.833682,
                 2.030204, 2.247787, 2.488689), tolerance = 1e-6)
})

context("3. test - locLinear2d: sparse case")
xcovResList <- lapply(rawCovList, function(x){
  xout <- seq(min(x$t1), max(x$t1), length.out = 30)
  x %>>% `[`(t1 != t2) %>>% {
    lapply(seq(10, 20, by = 5)/10, function(bw){
      locLinear2d(c(bw, bw), as.matrix(.[ , .(t1, t2)]), .$sse, .$weight, .$cnt, xout, xout, "gauss")
    })
  }
})
test_that("test - locLinear2d", {
  expect_equal(xcovResList[[1]][[1]], xcov_10_case2, tolerance = 1e-6)
  expect_equal(xcovResList[[1]][[2]], xcov_15_case2, tolerance = 1e-6)
  expect_equal(xcovResList[[1]][[3]], xcov_20_case2, tolerance = 1e-6)
  expect_equal(xcovResList[[2]][[1]], xcov_10_case1, tolerance = 1e-6)
  expect_equal(xcovResList[[2]][[2]], xcov_15_case1, tolerance = 1e-6)
  expect_equal(xcovResList[[2]][[3]], xcov_20_case1, tolerance = 1e-6)
  expect_equal(xcovResList[[3]][[1]], xcov_10_case0, tolerance = 1e-6)
  expect_equal(xcovResList[[3]][[2]], xcov_15_case0, tolerance = 1e-6)
  expect_equal(xcovResList[[3]][[3]], xcov_20_case0, tolerance = 1e-6)
})


