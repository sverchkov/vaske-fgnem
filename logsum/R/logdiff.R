"logdiff" <-
function(a, b, thresh=2e-13) {
  if (a < b) {
    warning("Log of a negative difference is NA")
    return(NaN)
  }
  if (a == b) {
    return(-Inf)
  }
  est1 <- b + log(expm1(a-b))
  est2 <- a + log(1 - exp(b-a))
 if (!is.finite(est1)) {
    warning("Bad estimate 1 ", as.character(est1),
            " with ", as.character(a), " and ", as.character(b))
    return(est2)
  }
   if (!is.finite(est2)) {
    warning("Bad estimate 2 ", as.character(est2),
            " with ", as.character(a), " and ", as.character(b))
    return(est1)
  }
  if (abs(est1 - est2) > thresh) {
    warning("Estimates of difference differ by more than threshold")
  }
  return (est1)
}

