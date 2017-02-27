\name{histratio}
\alias{histratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
histratio(x, breaks = "Sturges", conf.level = 0.95, length = 0.1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{breaks}{ ~~Describe \code{breaks} here~~ }
  \item{conf.level}{ ~~Describe \code{conf.level} here~~ }
  \item{length}{ ~~Describe \code{length} here~~ }
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
function(x, breaks="Sturges", conf.level=0.95, length=0.1, ...) {
  allhist <- hist(unlist(x), breaks = breaks, plot = FALSE)
  combhist <- sapply(x, function(z) hist(z, breaks = allhist$breaks,
                                           plot = FALSE)$counts)
  b <- binom.confint(combhist[,1], rowSums(combhist), methods="exact")
  x <- barplot(b[,'mean'], names=signif(allhist$mids, 2), ...)
  arrows(x, b[,'lower'], x, b[,'upper'], angle=90, code=3, length=length)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
