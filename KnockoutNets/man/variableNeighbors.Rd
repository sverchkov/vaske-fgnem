\name{variableNeighbors}
\alias{variableNeighbors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
variableNeighbors(scoreDF, edgeAssignments)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scoreDF}{ ~~Describe \code{scoreDF} here~~ }
  \item{edgeAssignments}{ ~~Describe \code{edgeAssignments} here~~ }
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
function (scoreDF, edgeAssignments) 
{
    stringIds <- idsOfScoreMatrix(scoreDF)
    ids <- matrix(match(as.matrix(scoreDF)[, 1:2], stringIds), 
        ncol = 2)
    result <- apply(ids, 1, function(x) {
        which(edgeAssignments[, 1] == x[1] & edgeAssignments[, 
            2] == x[2])
    })
    if (is.null(dim(result))) {
        dim(result) <- c(1, length(result))
    }
    return(result)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
