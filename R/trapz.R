#' trapz: trapezoidal rule to approximate the integral values
#'
#' Returns approximation of integral.
#'
#' @param x A vector with n elements, \code{x[i]} is a support, \code{i = 1, ..., n}.
#' If y is NULL, support is taken as \code{seq(1, length(x), by = 1)}.
#' @param y \code{y[i, j]} is jth values on corresponding value of \code{x[i], i = 1, ..., n}.
#' If y is vector, the length of y must be equal to the lenght of x.
#' If y is matrix, the number of rows must be equal to the lenght of x.
#' @return A value, the approximation of integral.
#' @section Reference:
#' Kai Habel, trapz, Octave.
#' @examples
#' # case 1
#' x <- c(1, 4, 9, 16, 25)
#' trapz(x) # 42
#'
#' # case 2
#' x <- matrix(c(1,4,9,16,25,1,8,27,64,125), 5)
#' trapz(x) # 42 162
#' @export
trapz <- function(x, y = NULL){
  if (is.null(y))
  {
    y <- x
    if (is.vector(x))
      x <- 1:length(x)
    else
      x <- 1:nrow(x)
  }
  return(trapz_cpp(x, as.matrix(y)))
}
