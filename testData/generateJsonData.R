set.seed(100)
library(pipeR)
library(rfda)
tp <- seq(1, 10, len = 101)
DT <- funcDataGen(100, tp, function(x) sin(x), function(x) rep(1, length(x)), "BesselJ")
sparseDT <- sparsify(DT, "sampleID", 0.85)
sparseDT %>>% `[`(j = .(t = list(t), y = list(y)), by = sampleID) %>>%
  toJSON %>>% write("../inst/extdata/funcData.json")
