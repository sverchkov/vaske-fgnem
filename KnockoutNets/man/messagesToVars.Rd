\name{messagesToVars}
\alias{messagesToVars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
messagesToVars(potential, msgsToFactors, norm = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{potential}{ ~~Describe \code{potential} here~~ }
  \item{msgsToFactors}{ ~~Describe \code{msgsToFactors} here~~ }
  \item{norm}{ ~~Describe \code{norm} here~~ }
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
function (potential, msgsToFactors, norm = T) 
{
    if (nrow(msgsToFactors)\%\%3 != 0 || ncol(msgsToFactors) != 
        6) {
        stop("malformed msgsToFactors")
    }
    result <- matrix(NA, nrow = nrow(msgsToFactors), ncol = ncol(msgsToFactors))
    for (i in seq(1, nrow(msgsToFactors), by = 3)) {
        ab <- msgsToFactors[i, ]
        ac <- msgsToFactors[i + 1, ]
        bc <- msgsToFactors[i + 2, ]
        bcOac <- outer(bc, ac, FUN = "+")
        result[i, ] <- apply(potential, 3, function(x) {
            max(x + bcOac)
        })
        bcOab <- outer(bc, ab, FUN = "+")
        result[i + 1, ] <- apply(potential, 2, function(x) {
            max(x + bcOab)
        })
        acOab <- outer(ac, ab, FUN = "+")
        result[i + 2, ] <- apply(potential, 1, function(x) {
            max(x + acOab)
        })
    }
    if (norm) {
        weights <- apply(result, 1, logsum)
    }
    else {
        weights <- apply(result, 1, median)
    }
    result <- result - weights
    return(result)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
