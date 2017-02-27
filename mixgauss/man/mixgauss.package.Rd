\name{mixgauss-package}
\alias{mixgauss-package}
\alias{mixgauss}
\docType{package}
\title{
Mixture of Gaussians
}
\description{
Basic mixture of Gaussians algorithm.  Given a set of one 
  dimensional points and a number of Gaussians, use EM to determine each 
  Gaussian's mean and standard deviation.
}
\details{
\tabular{ll}{
Package: \tab mixgauss\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-07-14\cr
License: \tab BSD\cr
}
Basic mixture of Gaussians algorithm.  Given a set of one 
  dimensional points and a number of Gaussians, use EM to determine each 
  Gaussian's mean and standard deviation.

}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
~~ Optionally other keywords from doc/KEYWORDS, one per line
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
  x <- c(rnorm(300, mean=0, sd=2), rnorm(500, mean=3, sd=5))
  mg <- mixGaussians(x, models=2, iter.max=15,
                     mu=c(-1,1), sigma=c(2,2), fixed.mu=c(0, NA))
}
