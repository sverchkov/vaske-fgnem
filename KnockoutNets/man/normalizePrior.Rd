\name{normalizePrior}
\alias{normalizePrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normalizePrior(prior, nullnull = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{ ~~Describe \code{prior} here~~ }
  \item{nullnull}{ ~~Describe \code{nullnull} here~~ }
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
function (prior, nullnull = 0.8) 
{
    prior[1] <- 1/3
    prior[2:length(prior)] <- (2/3)/(length(prior - 1) - 1)
    return(prior)
    return(prior/sum(prior))
    prior[1] <- nullnull
    norm <- (1 - nullnull)/sum(prior[2:length(prior)])
    prior[2:length(prior)] <- norm * prior[2:length(prior)]
    return(prior)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
