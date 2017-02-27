\name{edgeAssignments}
\alias{edgeAssignments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
edgeAssignments(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ ~~Describe \code{n} here~~ }
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
function (n) 
{
    result <- matrix(NA, nrow = (n * (n - 1) * (n - 2)/2), ncol = 3)
    idx <- 1
    for (i in 1:(n - 2)) {
        for (j in (i + 1):(n - 1)) {
            for (k in (j + 1):n) {
                if (i != k && j != k) {
                  result[idx + 0, ] <- c(i, j, k)
                  result[idx + 1, ] <- c(i, k, j)
                  result[idx + 2, ] <- c(j, k, i)
                  idx <- idx + 3
                }
            }
        }
    }
    return(result)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
