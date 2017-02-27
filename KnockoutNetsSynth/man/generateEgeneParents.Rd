\name{generateEgeneParents}
\alias{generateEgeneParents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generateEgeneParents(nEgenes, sgeneNames, prior = NULL, inhibition = 0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nEgenes}{ ~~Describe \code{nEgenes} here~~ }
  \item{sgeneNames}{ ~~Describe \code{sgeneNames} here~~ }
  \item{prior}{ ~~Describe \code{prior} here~~ }
  \item{inhibition}{ ~~Describe \code{inhibition} here~~ }
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
function (nEgenes, sgeneNames, prior=NULL,
                                  inhibition=0.25) {
  parent <- sample(1:length(sgeneNames), size=nEgenes,
                   replace=TRUE, prob=prior)
  sign <- sample(c(-1,1), size=nEgenes, replace=TRUE,
                 prob=c(inhibition, 1 - inhibition))
  parent <- sign * parent
  parentChar <- c("n", "error", "p")[sign+2]
  names(parent) <- paste('egene', 1:nEgenes, parentChar,
                         sgeneNames[abs(parent)], sep="")
  return(parent)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
