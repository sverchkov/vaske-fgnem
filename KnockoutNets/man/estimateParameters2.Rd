\name{estimateParameters2}
\alias{estimateParameters2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimateParameters2(egenes, EMiter = 15, mu = NULL, sigma = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egenes}{ ~~Describe \code{egenes} here~~ }
  \item{EMiter}{ ~~Describe \code{EMiter} here~~ }
  \item{mu}{ ~~Describe \code{mu} here~~ }
  \item{sigma}{ ~~Describe \code{sigma} here~~ }
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
function (egenes, EMiter = 15, mu = NULL, sigma = NULL) 
{
    x <- as.vector(egenes)
    x <- x[!is.na(x)]
    source("~/src/cvaske/R/MixGauss.R")
    mg <- mixGaussians(x, models = 2, iter.max = EMiter, mu = mu, 
        sigma = sigma, fixed.mu = c(0, NA))
    params <- cbind(mg$mu, mg$sigma, mg$alpha)
    colnames(params) <- c("mean", "sd", "alpha")
    if (sum(params[, "mean"]) > 0) {
        params <- rbind(c(Inf, 0, 0), params)
    }
    else {
        params <- rbind(c(Inf, 0, 0), params)[3:1, ]
    }
    rownames(params) <- c("neg", "null", "pos")
    return(params)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
