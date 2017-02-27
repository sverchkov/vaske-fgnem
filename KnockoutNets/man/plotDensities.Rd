\name{plotDensities}
\alias{plotDensities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotDensities(densa, densb, acol = 1, bcol = 2, labels = character(0), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{densa}{ ~~Describe \code{densa} here~~ }
  \item{densb}{ ~~Describe \code{densb} here~~ }
  \item{acol}{ ~~Describe \code{acol} here~~ }
  \item{bcol}{ ~~Describe \code{bcol} here~~ }
  \item{labels}{ ~~Describe \code{labels} here~~ }
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
function (densa, densb, acol = 1, bcol = 2, labels = character(0), 
    ...) 
{
    xlim = range(c(densa$x, densb$x))
    ylim = range(c(densa$y, densb$y))
    plot(densa, xlim = xlim, ylim = ylim, col = acol, ...)
    lines(densb$x, densb$y, col = bcol)
    if (length(labels) == 2) {
        legend(x = xlim[1], y = ylim[2], legend = labels, col = c(acol, 
            bcol), lty = 1)
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
