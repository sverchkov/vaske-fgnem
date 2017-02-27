#' Print Mixed Gaussians
#' @export
"print.mixGaussians" <-
function(x, ...) {
  cat("Mixture of ", length(x$mu), " Gaussians of ",
      nrow(x$y), " points.\n", sep="")
  cat("\nGeneration probabilities\n")
  l <- min(10, nrow(x$y))
  print(x$y[1:l,], ...)
  if (l < nrow(x$y)) {
    cat('... cutoff\n')
  }
  cat("\nDistribution means:\n")
  print(x$mu, ...)
  cat("\nDistribution standard deviations:\n")
  print(x$sigma, ...)
  cat("\nMixture Parameters:\n")
  print(x$alpha, ...)
}

