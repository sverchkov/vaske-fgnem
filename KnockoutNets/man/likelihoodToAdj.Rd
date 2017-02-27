\name{likelihoodToAdj}
\alias{likelihoodToAdj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
likelihoodToAdj(ll.df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ll.df}{ ~~Describe \code{ll.df} here~~ }
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
function (ll.df) 
{
    names <- unique(as.vector(t(ll.df[, 1:2])))
    r <- diag(length(names))
    dimnames(r) <- list(names, names)
    for (i in 1:nrow(ll.df)) {
        edge <- which.max(ll.df[i, 3:8])
        A <- ll.df[i, 1]
        B <- ll.df[i, 2]
        if (edge == 1) {
            r[A, B] <- 1
        }
        else if (edge == 2) {
            r[B, A] <- 1
        }
        else if (edge == 3) {
            r[A, B] <- -1
        }
        else if (edge == 4) {
            r[B, A] <- -1
        }
        else if (edge == 5) {
        }
        else if (edge == 6) {
            r[A, B] <- r[B, A] <- 1
        }
    }
    return(r)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
