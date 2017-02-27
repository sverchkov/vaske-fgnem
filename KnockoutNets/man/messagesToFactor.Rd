\name{messagesToFactor}
\alias{messagesToFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
messagesToFactor(constants, msgsToVars, varNeighbors, norm = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{constants}{ ~~Describe \code{constants} here~~ }
  \item{msgsToVars}{ ~~Describe \code{msgsToVars} here~~ }
  \item{varNeighbors}{ ~~Describe \code{varNeighbors} here~~ }
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
function (constants, msgsToVars, varNeighbors, norm = NULL) 
{
    if (nrow(msgsToVars) == 0 || ncol(msgsToVars) != ncol(constants)) {
        stop("malformed msgsToVars")
    }
    if (ncol(varNeighbors) != nrow(constants)) {
        stop("malformed varNeighbors")
    }
    result <- matrix(NA, nrow = nrow(msgsToVars), ncol = ncol(msgsToVars))
    for (i in 1:nrow(constants)) {
        outgoing <- varNeighbors[, i, drop = FALSE]
        total <- constants[i, ] + colSums(msgsToVars[outgoing, 
            , drop = FALSE])
        for (msg in outgoing) {
            result[msg, ] <- total - msgsToVars[msg, , drop = FALSE]
        }
    }
    if (!is.null(norm)) {
        result <- result - apply(result, 1, norm)
    }
    return(result)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
