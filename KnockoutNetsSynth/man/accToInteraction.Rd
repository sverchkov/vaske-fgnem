\name{accToInteraction}
\alias{accToInteraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
accToInteraction(acc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acc}{ ~~Describe \code{acc} here~~ }
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
function(acc) {
  kd <- colnames(acc)
  pairs <- expand.grid(A=kd,B=kd)[as.vector(upper.tri(diag(length(kd)))),]
  pair2int <- function(p) {
    ab <- acc[p[1],p[2]]
    ba <- acc[p[2],p[1]]
    i <- (sign(ab)+1) + 3*(sign(ba)+1) + 1
    map <- c(NA, # 1 - (-1, -1) - AeqvB - 6 NA
             4,  # 2 - (0,  -1) - BrepA - 4
             NA, # 3 - (1,  -1) - AeqvB - 6 NA
             3,  # 4 - (-1,  0) - ArepB - 3
             5,  # 5 - (0,   0) - AnonB - 5
             1,  # 6 - (1,   0) - AtoB  - 1
             NA, # 7 - (-1,  1) - AeqvB - 6 NA
             2,  # 8 - (0,   1) - BtoA  - 2
             6)  # 9 - (1,   1) - AeqvB - 6
#    map <- c(NA,4,NA,3,5,1,NA,2,6)
    map[i]
#    c(max(abs(ab), abs(ba)), map[i])
  }
  r <- apply(pairs, 1, pair2int)
#  rownames(r) <- apply(pairs, 1, paste, collapse="*")
  names(r) <- apply(pairs, 1, paste, collapse="*")
  r
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
