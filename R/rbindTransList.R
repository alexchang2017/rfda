#' rbind a transposed list of lists of data.frames
#'
#' @param x A list of lists of data.frames.
#' @return A list of data.tables.
#' @examples
#' ll <- list(list(data.frame(a = 1:2, b = 3:4), data.frame(a = 1:3, b = 3:5)),
#'            list(data.frame(a = 4:6, b = 2:4), data.frame(a = 3:2, b = 6:7)))
#' rbindTransList(ll)
#' @importFrom data.table rbindlist
#' @importFrom plyr laply ldply
#' @export
rbindTransList <- function(x){
  if (!is.list(x) || any(!laply(x, is.list)) || !all(laply(x, function(l) all(laply(l, is.data.frame)))))
    stop("x is not the list of lists of data.frames.")
  if (any(laply(x, length) != length(x[[1]])))
    stop("The lengths of elements of x is not all the same.")

  y <- vector('list', length(x[[1]]))
  for (i in 1:length(x[[1]]))
    y[[i]] <- ldply(x, `[[`, i) %>>% setDT
  return(y)
}
