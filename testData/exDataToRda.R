library(rmatio)
library(pipeR)
library(stringr)

# exData
matFiles <- list.files("exData", "\\.mat", full.names = TRUE)
exData <- matFiles %>>% lapply(read.mat) %>>%
  lapply(function(x){
    sampleID <- seq_along(x[[1]][[1]]) %>>%
      lapply(function(i) rep(i, length(x[[1]][[1]][[i]]))) %>>%
      do.call(what = c)
    data.frame(y = unlist(x[[1]][[1]]), t = unlist(x[[2]][[1]]), sampleID = sampleID)
  })

for (i in seq_along(exData)){
  dataName <- str_match(matFiles[i], "/([a-zA-z]*)\\.mat")[2]
  assign(dataName, exData[[i]])
  save(list = dataName, file = paste0("exData/", dataName,".rda"), compress = "xz")
}

# test result
matFiles <- list.files("testMat", "\\.mat", full.names = TRUE)
matFiles %>>% lapply(function(x){
  matRes <- read.mat(x)
  attach(matRes)
  savingName <- str_replace(x, "\\.mat", "\\.rda") %>>% str_replace("testMat/", "../tests/testthat/")
  save(list = names(matRes) , file = savingName, compress = "xz")
  detach(matRes)
})

