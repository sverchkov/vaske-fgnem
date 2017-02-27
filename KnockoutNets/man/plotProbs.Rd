\name{plotProbs}
\alias{plotProbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotProbs(scoreStruct, scaleLiklihood = F, zlim = c(-5, 5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scoreStruct}{ ~~Describe \code{scoreStruct} here~~ }
  \item{scaleLiklihood}{ ~~Describe \code{scaleLiklihood} here~~ }
  \item{zlim}{ ~~Describe \code{zlim} here~~ }
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
function (scoreStruct, scaleLiklihood = F, zlim = c(-5, 5)) 
{
    m <- sapply(scoreStruct$egenes.probs, function(x) x[, 3] - 
        x[, 1])
    if (scaleLiklihood) {
        r <- range(m)
        m[m < 0] <- m[m < 0]/-r[1]
        m[m > 0] <- m[m > 0]/r[2]
    }
    if (!is.null(zlim)) {
        m[m < zlim[1]] <- zlim[1]
        m[m > zlim[2]] <- zlim[2]
    }
    source("~/src/cvaske/R/plotting.R")
    yb.pal <- c(rgb(0, 0, 16:0/16), rgb(1:16/16, 1:16/16, 0))
    my.image(m, zlim = zlim, col = yb.pal)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
