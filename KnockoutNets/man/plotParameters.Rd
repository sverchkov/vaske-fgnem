\name{plotParameters}
\alias{plotParameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotParameters(egenes, p, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egenes}{ ~~Describe \code{egenes} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
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
function (egenes, p, ...) 
{
    egenes <- egenes[!is.na(egenes)]
    mesh <- seq(from = range(egenes)[1], to = range(egenes)[2], 
        length.out = 500)
    hist(as.vector(egenes), prob = TRUE, main = "Egene observations", 
        ...)
    linecol <- c(3, 1, 2)[rank(p[, "mean"])]
    for (i in 1:nrow(p)) {
        lines(mesh, dnorm(mesh, mean = p[i, "mean"], sd = p[i, 
            "sd"]) * p[i, "alpha"], col = linecol[i])
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
