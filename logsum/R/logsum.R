#' @export
"logsum" <-
function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA)
  }
  if (length(x) < 2) {
    return(x)
  }
  result <- logsum.two(x[1], x[2])
  if (length(x) > 2) {
    for (i in 3:length(x)) {
      result <- logsum.two(result, x[i])
    }
  }
  return(result)
}

