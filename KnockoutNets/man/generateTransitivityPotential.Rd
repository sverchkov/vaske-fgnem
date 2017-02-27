\name{generateTransitivityPotential}
\alias{generateTransitivityPotential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generateTransitivityPotential()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    n <- c("->", "<-", "-|", "|-", ":", "<->")
    f <- c(1, 0, -1, 0, 0, 1)
    b <- c(0, 1, 0, -1, 0, 1)
    result <- matrix(rep(NA, 6 * 6 * 6))
    dim(result) <- c(6, 6, 6)
    dimnames(result) <- list(paste("B", n, "C", sep = ""), paste("A", 
        n, "C", sep = ""), paste("A", n, "B", sep = ""))
    for (i in 1:6) {
        Mab <- f[i]
        Mba <- b[i]
        for (j in 1:6) {
            Mac <- f[j]
            Mca <- b[j]
            for (k in 1:6) {
                Mbc <- f[k]
                Mcb <- b[k]
                notTrans <- function(ab, bc, ac) {
                  return(ab & bc & (ac != ab * bc))
                }
                chains <- matrix(c(Mab, Mbc, Mac, Mac, Mcb, Mab, 
                  Mba, Mac, Mbc, Mbc, Mca, Mba, Mca, Mab, Mcb, 
                  Mcb, Mba, Mca), nrow = 3)
                mismatch <- any(apply(chains, 2, function(x) notTrans(x[1], 
                  x[2], x[3])))
                result[i, j, k] <- if (mismatch) 
                  -Inf
                else 0
            }
        }
    }
    return(result)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
