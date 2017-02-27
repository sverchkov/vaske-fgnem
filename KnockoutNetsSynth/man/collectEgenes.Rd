\name{collectEgenes}
\alias{collectEgenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
collectEgenes(data, posterior, scorematrix = "scores")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{posterior}{ ~~Describe \code{posterior} here~~ }
  \item{scorematrix}{ ~~Describe \code{scorematrix} here~~ }
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
function(data, posterior, scorematrix="scores") {
  stopifnot(length(data)==length(posterior))
  scores <- lapply(1:length(data), function(i) {
    scorePairwise(predictions[[i]][[scorematrix]],
                  data[[i]], data[[i]]$obsParams)
  })
  negenes <- sapply(scores, function(x) nrow(x$posterior))
  id <- rep(1:length(data), negenes)
  fracInhib <- rep(sapply(data, function(d) sum(d$sgeneAdj < 0) /
                          sum(d$sgeneAdj != 0)), negenes)
  ll <- rep(sapply(scores, function(x) x$ll), negenes)
  actualEgeneParent <- unlist(lapply(data, function(d) d$egeneParents))
  predEgeneParent <- unlist(lapply(scores, egenePostToParent))
  cbind(id, fracInhib, ll, actualEgeneParent, predEgeneParent)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
