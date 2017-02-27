\name{generateExpression}
\alias{generateExpression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generateExpression(quant, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{quant}{ ~~Describe \code{quant} here~~ }
  \item{params}{ ~~Describe \code{params} here~~ }
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
function(quant, params) {
  r <- matrix(NA, nrow=nrow(quant), ncol=ncol(quant), dimnames=dimnames(quant))
  r[quant == -1] <- rnorm(sum(quant == -1), mean=params["neg", "mean"],
                          sd=params["neg", "sd"])
  r[quant == 0] <- rnorm(sum(quant == 0), mean=params["null", "mean"],
      sd=params["null", "sd"])
  r[quant == 1] <- rnorm(sum(quant == 1), mean=params["pos", "mean"],
      sd=params["pos", "sd"])
  return(r)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
