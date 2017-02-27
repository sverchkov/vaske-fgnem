\name{logsum}
\alias{logsum}
\alias{logsum.two}
\title{ logsum, logsum.two }
\description{
  Add numbers that are in logspace.
}
\usage{
logsum(x)
logsum.two(a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector of numbers }
  \item{a}{ a single number }
  \item{b}{ a single number }

}
\details{
  Adds two numbers in logspace
}
\value{
  The addition of it's arguments, also in logspace
}
\references{ Durbin, Eddy, Krogh, Mitchison, Biological Sequence
  Analysis, Section 3.6 }
\author{ Charlie Vaske }
\seealso{ See Also as \code{\link{logdiff}} }
\examples{
logsum(c(-Inf, -3))
logsum.two(-Inf, -3)
logsum(-10:-1)

d <- logsum(c(-3, -5))
e <- log(exp(-3) + exp(-5))
stopifnot(d - e < 1e-10)
}
\keyword{ arith }
