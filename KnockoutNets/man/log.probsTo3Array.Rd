\name{log.probsTo3Array}
\alias{log.probsTo3Array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
log.probsTo3Array(egenes.logprobs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egenes.logprobs}{ ~~Describe \code{egenes.logprobs} here~~ }
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
function (egenes.logprobs) 
{
    array(unlist(egenes.logprobs), c(dim(egenes.logprobs[[1]]), 
        length(egenes.logprobs)), dimnames = list(rownames(egenes.logprobs[[1]]), 
        colnames(egenes.logprobs[[1]]), names(egenes.logprobs)))
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
