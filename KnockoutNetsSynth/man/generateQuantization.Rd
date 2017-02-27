\name{generateQuantization}
\alias{generateQuantization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generateQuantization(knockdown.cols, parents, acc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{knockdown.cols}{ ~~Describe \code{knockdown.cols} here~~ }
  \item{parents}{ ~~Describe \code{parents} here~~ }
  \item{acc}{ ~~Describe \code{acc} here~~ }
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
function(knockdown.cols, parents, acc) {
  sapply(knockdown.cols, function(x) {
    p <- acc[x,abs(parents)] * sign(parents)
    r <- (runif(parents) < abs(p)) * sign(p)
    names(r) <- names(parents)
    r
  })
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
