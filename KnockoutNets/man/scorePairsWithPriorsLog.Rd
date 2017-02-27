\name{scorePairsWithPriorsLog}
\alias{scorePairsWithPriorsLog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scorePairsWithPriorsLog(pairs, egenes.logprobs, params, summarization)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pairs}{ ~~Describe \code{pairs} here~~ }
  \item{egenes.logprobs}{ ~~Describe \code{egenes.logprobs} here~~ }
  \item{params}{ ~~Describe \code{params} here~~ }
  \item{summarization}{ ~~Describe \code{summarization} here~~ }
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
function (pairs, egenes.logprobs, params, summarization) 
{
    AtoB <- apply(pairs, 1, function(x) {
        sGenePairDomLog(egenes.logprobs[[x[1]]], egenes.logprobs[[x[2]]], 
            params, summarization = summarization)
    })
    BtoA <- apply(pairs, 1, function(x) {
        sGenePairDomLog(egenes.logprobs[[x[2]]], egenes.logprobs[[x[1]]], 
            params, summarization = summarization)
    })
    ArepB <- apply(pairs, 1, function(x) {
        sGenePairRepLog(egenes.logprobs[[x[1]]], egenes.logprobs[[x[2]]], 
            params, summarization = summarization)
    })
    BrepA <- apply(pairs, 1, function(x) {
        sGenePairRepLog(egenes.logprobs[[x[2]]], egenes.logprobs[[x[1]]], 
            params, summarization = summarization)
    })
    AnonB <- apply(pairs, 1, function(x) {
        sGenePairNonLog(egenes.logprobs[[x[1]]], egenes.logprobs[[x[2]]], 
            params, summarization = summarization)
    })
    AeqvB <- apply(pairs, 1, function(x) {
        sGenePairEqvLog(egenes.logprobs[[x[1]]], egenes.logprobs[[x[2]]], 
            params, summarization = summarization)
    })
    ll <- cbind(pairs, AtoB = AtoB, BtoA = BtoA, ArepB, BrepA, 
        AnonB, AeqvB)
    ll
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
