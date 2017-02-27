\name{logdiff}
\alias{logdiff}
\title{ ~~function to do ... ~~ }
\description{
  Subtract \code{b} from \code{a}, for numbers in log space
}
\usage{
logdiff(a, b, thresh = 2e-13)
}
\arguments{
  \item{a}{ the larger number }
  \item{b}{ the smaller number }
  \item{thresh}{ threshold for reporting inaccurate estimates }
}
\details{
  Subtracts \code{a} from \code{b}, for numbers that are in logspace.
  Two different estimates are made, and compared for accuracy.
}
\value{
  The difference of \code{a} and \code{b}, in log space.
}
\author{ Charlie Vaske }
\seealso{ See also \code{\link{logsum}} }
\examples{
d <- logdiff(-3, -5)

e <- log(exp(-3) - exp(-5))

stopifnot(d - e < 1e-10)
}
\keyword{ arith }% at least one, from doc/KEYWORDS

