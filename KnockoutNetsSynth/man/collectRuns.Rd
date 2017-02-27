\name{collectRuns}
\alias{collectRuns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
collectRuns(data, predictions, scorematrix = "scores")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{predictions}{ ~~Describe \code{predictions} here~~ }
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
function(data, predictions, scorematrix="scores") {
  if (length(data) != length(predictions)) {
    stop("length mismatch on data and prediction lists")
  }
  tmp.sizes <- sapply(data, function(x) nrow(x$sgeneAdj))
  id <- rep(1:length(data), (tmp.sizes -1) * tmp.sizes / 2)
  size <- unlist(lapply(data, gatherSize))
  separation <- unlist(lapply(data, gatherMean))
  fracInhib <- unlist(lapply(data, gatherFracInhibition))
  strength <- unlist(lapply(data, function(d) accToLinkStrength(d$sgeneAcc)))
  actual <- unlist(lapply(data, function(d) accToInteraction(d$sgeneAcc)))
  predicted <- unlist(lapply(predictions, function(p)
                             likelihoodToPredictions(p[[scorematrix]])))
  cbind(id, size, separation, fracInhib, strength, actual, predicted)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
