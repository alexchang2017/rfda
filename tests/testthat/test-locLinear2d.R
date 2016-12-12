context("test - locLinear2d")

load("covRes.rda") # results from MatLab
data("regularExData", package = "rfda")
data("irregularExData", package = "rfda")
data("sparseExData", package = "rfda")

testDataRegular <- data.table(y = c(seq(0.1, 0.9, 0.1), c(0.9,1,1)), t = rep(1:3, 4),
                              sampleID = c(rep(1:4, each=3)), key = "t")
testDataIrregular <- data.table(y = seq(0.1, 1.3, 0.1), t = c(1:4,1,2,4,2:4,1,3,4),
                                sampleID = c(rep(1, 4), rep(2:4, each=3)), key = "t")
testDataSparse <- data.table(y = seq(0.1, 1, 0.1), t = c(1,2,4,3,5,7,6,10,8,9),
                             sampleID = c(rep(1:2, each=3), rep(3:4, each=2)), key = "t")
testDT_list <- suppressMessages(lapply(list(testDataRegular, testDataIrregular, testDataSparse), function(dt){
  sparsity <- checkSparsity(dt, "sampleID", "t")
  bwCand <- bwCandChooser(dt, "sampleID", "t", sparsity, "gauss", 1)
  w <- rep(1, nrow(dt))
  bwOpt <- gcvLocPoly1d(bwCand, dt$t, dt$y, w, "gauss", 0, 1)
  bwOpt <- adjGcvBw(bwOpt, sparsity, "gauss", 0)
  xout <- sort(unique(dt$t))
  meanFunc <- locPoly1d(bwOpt, dt$t, dt$y, w, xout, "gauss", 0, 1)

  merge(dt, data.table(mf = meanFunc, t = xout, key = "t"), by = "t") %>>%
    `[`(j = `:=`(y = y - mf, variable = "y")) %>>%
    setnames(c("t", "y", "sampleID"), c("timePnt", "value.demean", "subId")) %>>%
    `[`(j = mf := NULL)
}))

compareCovResList <- lapply(list(rcov_case2, rcov_case1, rcov_case0), function(m){
  data.table(m) %>>% setnames(c("t1", "t2", "sse")) %>>% setorder(t1, t2)
})

context("1. test - getCrCov")
rawCovList <- lapply(testDT_list, getRawCrCov)
test_that("test - bwCandChooser2", {
  expect_equal(rawCovList[[1]][ , .(t1,t2,sse)], compareCovResList[[1]] %>>% `[`(j = sse := sse * 4))
  expect_equal(rawCovList[[2]][ , .(t1,t2,sse)], compareCovResList[[2]])
  expect_equal(rawCovList[[3]][ , .(t1,t2,sse)], compareCovResList[[3]])
})

context("2. test - bwCandChooser2")
allGridList <- lapply(list(regularExData, irregularExData, sparseExData), function(df){
  data.table(df) %>>% `[`(j = .(t1 = rep(t, length(t)), t2 = rep(t, each=length(t))), by = .(sampleID))
})
testErrDt <- data.table(sampleID = rep(1:4, each=2), t = c(1,3,11,9,1,13,4,2), y = rnorm(8)) %>>%
  `[`(j = .(t1 = rep(t, length(t)), t2 = rep(t, each=length(t))), by = .(sampleID))
test_that("test - bwCandChooser2", {
  expect_equal(bwCandChooser2(allGridList[[1]], 2, "gauss", 1)[ , 1],
               c(0.315789, 0.397407, 0.500119, 0.629378, 0.792045, 0.996754, 1.254371,
                 1.578570, 1.986561, 2.500000), tolerance = 1e-6)
  expect_equal(bwCandChooser2(allGridList[[2]], 1, "gauss", 1)[ , 1],
               c(0.210526, 0.277147, 0.364850, 0.480306, 0.632297, 0.832387, 1.095794,
                 1.442556, 1.899049, 2.500000), tolerance = 1e-6)
  expect_equal(bwCandChooser2(allGridList[[3]], 0, "gauss", 1)[ , 1],
               c(0.995476, 1.102164, 1.220286, 1.351068, 1.495867, 1.656184, 1.833682,
                 2.030204, 2.247787, 2.488689), tolerance = 1e-6)
  expect_message(bwCandChooser2(testErrDt, 0, "gauss", 9), "Data is too sparse")
  expect_error(bwCandChooser2(testErrDt, 0, "epan", 9), "Data is too sparse")
})

context("3. test - locLinear2d")
xcovResList <- lapply(rawCovList, function(x){
  xout <- seq(min(x$t1), max(x$t1), length.out = 30)
  x %>>% `[`(t1 != t2) %>>% {
    lapply(seq(10, 20, by = 5)/10, function(bw){
      locLinear2d(c(bw, bw), as.matrix(.[ , .(t1, t2)]), .$sse, .$weight, .$cnt, xout, xout, "gauss")
    })
  }
})
rcov_test <- rawCovList[[1]]
xout <- seq(min(rcov_test$t1), max(rcov_test$t1), length.out = 30)
xcovResList2 <- lapply(c("epan", "gaussvar", "quar"), function(k){
  rcov_test %>>% `[`(t1 != t2) %>>% {
    locLinear2d(c(3, 3), as.matrix(.[ , .(t1, t2)]), .$sse, .$weight, .$cnt, xout, xout, k)
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
  expect_equal(xcovResList2[[1]], xcov_30_epan, tolerance = 1e-6)
  expect_equal(xcovResList2[[3]], xcov_30_quar, tolerance = 1e-6)
  expect_equal(xcovResList2[[2]], xcov_30_gaussvar, tolerance = 1e-5) # there is a percision issue
})

context("4. test - locLinear2d: validate inputs")
m <- rcov_test[t1 != t2]
test_that("test - locLinear2d: validate inputs", {
  expect_true(is.na(locLinear2d(c(0.3, 0.3), as.matrix(m[ , .(t1, t2)]), m$sse,
                                m$weight, m$cnt, xout, xout, "epan")))
  expect_error(locLinear2d(c(3, -2), as.matrix(m[ , .(t1, t2)]), m$sse, m$weight, m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, NA), as.matrix(m[ , .(t1, t2)]), m$sse, m$weight, m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, NaN), as.matrix(m[ , .(t1, t2)]), m$sse, m$weight, m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, Inf), as.matrix(m[ , .(t1, t2)]), m$sse, m$weight, m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3, 3), as.matrix(m[ , .(t1, t2)]), m$sse, m$weight, m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2, cnt)]), m$sse,
                           m$weight, m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), m$sse,
                           m$weight, m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, NA), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(1, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, NaN), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(1, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, Inf), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(1, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), c(1, m$sse),
                           m$weight, m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(NA, m$sse),
                           c(1, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(NaN, m$sse),
                           c(1, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(Inf, m$sse),
                           c(1, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           c(1, m$weight), m$cnt, xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(NA, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(NaN, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(Inf, m$weight), c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, c(1, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(1, m$weight), c(NA, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(1, m$weight), c(NaN, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), rbind(c(1, 1), as.matrix(m[ , .(t1, t2)])), c(1, m$sse),
                           c(1, m$weight), c(Inf, m$cnt), xout, xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, c(2, xout), xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, c(NA, xout), xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, c(NaN, xout), xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, c(Inf, xout), xout, "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, xout, c(2, xout), "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, xout, c(NA, xout), "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, xout, c(NaN, xout), "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, xout, c(Inf, xout), "gauss"))
  expect_error(locLinear2d(c(3, 3), as.matrix(m[ , .(t1, t2)]), m$sse,
                           m$weight, m$cnt, xout, xout, "guass"))
})

context("5. test - gcvLocLinear2d check results")
sparsities <- lapply(testDT_list, checkSparsity, id.var = "subId", timeVarName = "timePnt")
gcvBws <- suppressMessages(mapply(function(rawCov, sparsity){
  bwCand <- bwCandChooser2(rawCov, sparsity, "gauss", 1)
  rawCovNoDiag <- rawCov[t1 != t2]
  gcvLocLinear2d(bwCand, as.matrix(rawCovNoDiag[ , .(t1, t2)]), rawCovNoDiag$sse,
                 rawCovNoDiag$weight, rawCovNoDiag$cnt, "gauss")
}, rawCovList, sparsities))
bwCand_test <- bwCandChooser2(rawCovList[[1]], 2, "epan", 1)
rawCovNoDiag_test <- rawCovList[[1]][t1 != t2]
test_that("test - gcvLocLinear2d check results", {
  expect_equal(gcvBws[ , 1], gcv_bw_case2, tolerance = 1e-6)
  expect_equal(gcvBws[ , 2], gcv_bw_case1, tolerance = 1e-6)
  expect_equal(gcvBws[ , 3], gcv_bw_case0, tolerance = 1e-6)
  expect_equal(with(rawCovNoDiag_test, gcvLocLinear2d(matrix(c(0.3, 0.5, 0.3, 0.5), 2), cbind(t1, t2),
                                                      sse, weight, cnt, "gaussvar")), c(0.3, 0.3), tolerance = 1e-6)
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse, weight, cnt, "epan")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse, weight, cnt, "quar")))
})

context("6. test - gcvLocLinear2d validate inputs")
testDT <- data.table(variable = rep(1, 6), timePnt = c(0, 10, 11, 100, 2, 5),
                     subId = rep(1:3, each = 2), value.demean = 1:6)
testRawCov <- getRawCrCov(testDT)
testBbwList <- bwCandChooser2(testRawCov, 0, "gauss", 1)
test_that("test - gcvLocLinear2d validate inputs", {
  expect_error(with(testRawCov, gcvLocLinear2d(testBbwList, cbind(t1, t2), sse, weight, cnt, "gauss")),
               "You may want to change to Gaussian kernel!")
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(rbind(c(-1, 1), bwCand_test), cbind(t1, t2), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(rbind(c(0, 1), bwCand_test), cbind(t1, t2), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(rbind(c(NA, 1), bwCand_test), cbind(t1, t2), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(rbind(c(NaN, 1), bwCand_test), cbind(t1, t2), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(rbind(c(Inf, 1), bwCand_test), cbind(t1, t2), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, rbind(c(NA, 1), cbind(t1, t2)), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, rbind(c(NaN, 1), cbind(t1, t2)), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, rbind(c(Inf, 1), cbind(t1, t2)), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, rbind(c(1, 1), cbind(t1, t2)), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2, t2), sse,
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), c(NA, sse),
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), c(NaN, sse),
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), c(Inf, sse),
                                                      weight, cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      c(NA, weight), cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      c(NaN, weight), cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      c(Inf, weight), cnt, "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, c(NA, cnt), "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, c(NaN, cnt), "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, c(Inf, cnt), "gauss")))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, cnt, "gauss", NA)))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, cnt, "gauss", NaN)))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, cnt, "gauss", Inf)))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, cnt, "gauss", -1)))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, cnt, "gauss", 0)))
  expect_error(with(rawCovNoDiag_test, gcvLocLinear2d(bwCand_test, cbind(t1, t2), sse,
                                                      weight, cnt, "gauss", 2.5)))
})

context("7. test - adjGcvBw")
test_that("test - adjGcvBw", {
  expect_equal(adjGcvBw(c(0.5, 0.5), 2, "gauss", 0), c(0.55, 0.55), tolerance = 1e-6)
  expect_equal(adjGcvBw(c(0.5, 0.5), 2, "epan", 0), c(0.55, 0.55), tolerance = 1e-6)
  expect_equal(adjGcvBw(c(0.4, 0.4), 1, "gauss", 0), c(0.44, 0.44), tolerance = 1e-6)
  expect_equal(adjGcvBw(c(0.4, 0.4), 1, "epan", 0), c(0.44, 0.44), tolerance = 1e-6)
  expect_equal(adjGcvBw(c(0.9964559, 0.9964559), 0, "gauss", 0), c(1.096101, 1.096101), tolerance = 1e-6)
  expect_equal(adjGcvBw(c(0.9964559, 0.9964559), 0, "epan", 0), c(1.096101, 1.096101), tolerance = 1e-6)
})

context("8. test - bwCandChooser3")
testDataRegularMulti <- data.table(y = c(seq(0.1, 0.9, 0.1), c(0.9,1,1)), t = rep(1:3, 4),
                                   sampleID = c(rep(1:4, each=3)), key = "t") %>>%
  `[`(j = `:=`(y2 = y * sin(t), y3 = y * cos(t)))
testDataIrregularMulti <- data.table(y = seq(0.1, 1.3, 0.1), t = c(1:4,1,2,4,2:4,1,3,4),
                                     sampleID = c(rep(1, 4), rep(2:4, each=3)), key = "t") %>>%
  `[`(j = `:=`(y2 = y * sin(t), y3 = y * cos(t)))
testDataSparseMulti <- data.table(y = seq(0.1, 1, 0.1), t = c(1,2,4,3,5,7,6,10,8,9),
                                  sampleID = c(rep(1:2, each=3), rep(3:4, each=2)), key = "t") %>>%
  `[`(j = `:=`(y2 = y * sin(t), y3 = y * cos(t)))

testDT_list2 <- suppressMessages(lapply(list(testDataRegularMulti, testDataIrregularMulti,
                                            testDataSparseMulti), function(dt){
  varNames <- c("y", "y2", "y3")
  sparsity <- checkSparsity(dt, "sampleID", "t")
  bwCand <- bwCandChooser(dt, "sampleID", "t", sparsity, "gauss", 1)
  w <- rep(1, nrow(dt))
  meanFuncList <- lapply(varNames, function(var){
    bwOpt <- gcvLocPoly1d(bwCand, dt$t, dt[[var]], w, "gauss", 0, 1)
    bwOpt <- adjGcvBw(bwOpt, sparsity, "gauss", 0)
    xout <- sort(unique(dt$t))
    data.table(t = xout, meanFunc = locPoly1d(bwOpt, dt$t, dt[[var]], w, xout, "gauss", 0, 1), key = "t") %>>%
      setnames("meanFunc", paste0("meanFunc_", var))
  })

  meanDT <- meanFuncList[[1]]
  for (i in 2:length(meanFuncList))
    meanDT <- merge(meanDT, meanFuncList[[i]], by = "t")

  expr <- paste0(varNames, "=", varNames, "-meanFunc_", varNames) %>>%
    paste0(collapse = ",") %>>% (paste0("`:=`(", ., ")"))
  merge(dt, meanDT, by = "t") %>>% `[`(j = eval(parse(text = expr))) %>>% `[`(j = variable := "y") %>>%
    melt.data.table(id.var = c("sampleID", "t"), measure.vars = varNames, variable.factor = FALSE, value.name = "value.demean") %>>%
    setnames(c("t", "sampleID"), c("timePnt", "subId"))
}))
rawCovList2 <- lapply(testDT_list2, getRawCrCov)
testErrDt <- data.table(sampleID = rep(1:4, each=2), t = c(1,3,11,9,1,13,4,2), y = rnorm(8)) %>>%
  `[`(j = .(t1 = rep(t, length(t)), t2 = rep(t, each=length(t))), by = .(sampleID))
test_that("test - bwCandChooser3", {
  expect_equal(unique(bwCandChooser3(rawCovList2[[1]], 2, "gauss", 1)[ , 1]), 0.5, tolerance = 1e-6)
  expect_equal(unique(bwCandChooser3(rawCovList2[[2]], 1, "gauss", 1)[ , 1]),
               c(0.4, 0.4680695, 0.5477226, 0.6409305, 0.75), tolerance = 1e-6)
  expect_equal(unique(bwCandChooser3(rawCovList2[[3]], 0, "gauss", 1)[ , 1]),
               c(0.9, 1.131690, 1.423025, 1.789359, 2.25), tolerance = 1e-6)
  expect_message(bwCandChooser3(testErrDt, 0, "gauss", 9), "Data is too sparse")
  expect_error(bwCandChooser3(testErrDt, 0, "epan", 9), "Data is too sparse")
})

context("9. test - locLinearRotate2d check results")
x <- CJ(seq(0, 1, 0.1), seq(0, 1, 0.1)) %>>% `[`(V1 != V2) %>>% as.matrix
idx <- x[,1] > x[,2]
f <- function(x) x ** 2 *2 + 3 * x + 3
y <- f(x[,1]) * (f(x[,2]) * 0.5 - 1)
w <- rep(1, length(y))
w2 <- rep(1:2, each = length(y) / 2)
count <- rep(1, length(y))
count2 <- rep(1:2, each = length(y) / 2)
outMat <- cbind(seq(0, 1, 0.2), seq(0, 1, 0.2))
test_that("test - locLinearRotate2d check results", {
  expect_equal(rfda:::locLinearRotate2d(c(1, 1), x, y, w, count, outMat, "gauss"),
               c(-1.125232, 2.734192, 6.745392, 10.910800, 15.233159, 19.715540), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(2, 2), x, y, w, count, outMat, "gauss"),
               c(-1.431382, 2.653423, 6.778417, 10.943776, 15.149684, 19.396329), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(1, 1), x, y, w, count, outMat, "gaussvar"),
               c(-0.798469, 2.804388, 6.710836, 10.875476, 15.303718, 20.053595), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(2, 2), x, y, w, count, outMat, "gaussvar"),
               c(-1.356083, 2.672359, 6.770155, 10.935285, 15.168108, 19.471385), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(1, 1), x, y, w, count, outMat, "epan"),
               c(0.247116, 2.962383, 6.675713, 10.845262, 15.455296, 21.215358), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(2, 2), x, y, w, count, outMat, "epan"),
               c(-1.274003, 2.687789, 6.766201, 10.931756, 15.183759, 19.549635), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(1, 1), x, y, w, count, outMat, "quar"),
               c(0.568979, 3.077283, 6.569519, 10.740611, 15.612548, 21.772818), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(2, 2), x, y, w, count, outMat, "quar"),
               c(-1.032552, 2.745381, 6.742974, 10.908719, 15.243003, 19.799011), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(1, 1), x, y, w, count2, outMat, "gauss"),
               c(1.369461, 3.021149, 4.679596, 6.351233, 8.043005, 9.762371), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(1, 1), x, y, w2, count, outMat, "gauss"),
               c(-1.516499, 2.466812, 6.613778, 10.929288, 15.418294, 20.085685), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(1, 1), x, y, w, count2, outMat, "epan"),
               c(1.371086, 3.091233, 4.673803, 6.314644, 8.002719, 9.938969), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(1, 1), x, y, w2, count, outMat, "epan"),
               c(0.229353, 2.833910, 6.571697, 10.847529, 15.587139, 21.489738), tolerance = 1e-6)
  expect_equal(rfda:::locLinearRotate2d(c(0.075, 0.075), x[idx, ], y[idx], w[idx], count[idx], outMat, "epan"),
               c(1.660000, 2.927927, 5.500262, 9.490163, 15.377616, 21.280000), tolerance = 1e-6)
})

context("10. test - locLinearRotate2d validate inputs")
test_that("test - locLinearRotate2d: validate inputs", {
  expect_true(is.na(rfda:::locLinearRotate2d(c(0.05, 0.05), x, y, w, count, outMat, "epan")))
  expect_error(rfda:::locLinearRotate2d(c(3, -2), x, y, w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, NA), x, y, w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, NaN), x, y, w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, Inf), x, y, w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3, 3), x, y, w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3, 3), x, y, w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), cbind(1, x), y, w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), y, w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, NA), x), c(1, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, NaN), x), c(1, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, Inf), x), c(1, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), x, c(1, y), w, count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(NA, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(NaN, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(Inf, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), x, y, c(1, w), count, outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(NA, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(NaN, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(Inf, y), c(1, 2), c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), x, y, w, c(1, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(1, y), c(1, 2), c(NA, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(1, y), c(1, 2), c(NaN, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), rbind(c(1, 1), x), c(1, y), c(1, 2), c(Inf, count), outMat, "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), x, y, w, count, cbind(1, outMat), "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), x, y, w, count, rbind(c(1, NA), outMat), "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), x, y, w, count, rbind(c(1, NaN), outMat), "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), x, y, w, count, rbind(c(1, Inf), outMat), "gauss"))
  expect_error(rfda:::locLinearRotate2d(c(3, 3), x, y, w, count, outMat, "guass"))
})
