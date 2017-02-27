\name{estimateParameters}
\alias{estimateParameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimateParameters(egenes, EMiter = 15, mu = NULL, sigma = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egenes}{ ~~Describe \code{egenes} here~~ }
  \item{EMiter}{ ~~Describe \code{EMiter} here~~ }
  \item{mu}{ ~~Describe \code{mu} here~~ }
  \item{sigma}{ ~~Describe \code{sigma} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
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
function (egenes, EMiter = 15, mu = NULL, sigma = NULL, verbose = FALSE) 
{
    source("~/src/cvaske/R/MixGauss.R")
    x <- as.vector(egenes)
    x <- x[!is.na(x)]
    mg <- mixGaussians(x, models = 3, iter.max = EMiter, mu = mu, 
        sigma = sigma, fixed.mu = c(NA, 0, NA))
    if (verbose) {
        print(mg)
    }
    params <- cbind(mg$mu, mg$sigma, mg$alpha)
    params <- params[order(mg$mu), ]
    colnames(params) <- c("mean", "sd", "alpha")
    rownames(params) <- c("neg", "null", "pos")
    return(params)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
