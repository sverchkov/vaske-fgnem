\name{print.mixGaussians}
\alias{print.mixGaussians}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.mixGaussians(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
