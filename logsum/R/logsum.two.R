"logsum.two" <-
function(a,b) {
  s <- a < b
  if (is.na(s)) {
    return(NA)
  } else if (a < b) {
    r <- c(a, b)
  } else {
    r <- c(b, a)
  }
  if (r[1] == -Inf) {
    return(r[2])
  } else {
    return(r[2] + log1p(exp(r[1]-r[2])))
  }
}

