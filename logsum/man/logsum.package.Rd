\name{logsum-package}
\alias{logsum-package}
\docType{package}
\title{
Basic arithmetic on numbers in log space.
}
\description{
Addition and subtraction for numbers that are in logspace
}
\details{
  Attempts to provide mostly numerically stable methods for addition and
  subtraction of log space numbers.
}
\author{
Charlie Vaske

Maintainer: Charlie Vaske <cvaske@soe.ucsc.edu>
}
\references{
Durbin, Eddy, Krogh, Mitchison, Biological Sequence Analysis, Section 3.6
}
\keyword{ package }
\examples{
logsum(c(-Inf, -3))
logsum.two(-Inf, -3)
logsum(-10:-1)
logdiff(-3, -5)
}
