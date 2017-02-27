\name{scorePairsWithPriors}
\alias{scorePairsWithPriors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scorePairsWithPriors(pairs, egenes.probs, egenes.logprobs, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pairs}{ ~~Describe \code{pairs} here~~ }
  \item{egenes.probs}{ ~~Describe \code{egenes.probs} here~~ }
  \item{egenes.logprobs}{ ~~Describe \code{egenes.logprobs} here~~ }
  \item{params}{ ~~Describe \code{params} here~~ }
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
function (pairs, egenes.probs, egenes.logprobs, params) 
{
    AtoB <- apply(pairs, 1, function(x) {
        sGenePairDom(egenes.probs[[x[1]]], egenes.probs[[x[2]]], 
            params)
    })
    BtoA <- apply(pairs, 1, function(x) {
        sGenePairDom(egenes.probs[[x[2]]], egenes.probs[[x[1]]], 
            params)
    })
    ArepB <- apply(pairs, 1, function(x) {
        sGenePairRep(egenes.probs[[x[1]]], egenes.probs[[x[2]]], 
            params)
    })
    BrepA <- apply(pairs, 1, function(x) {
        sGenePairRep(egenes.probs[[x[2]]], egenes.probs[[x[1]]], 
            params)
    })
    AnonB <- apply(pairs, 1, function(x) {
        sGenePairNon(egenes.probs[[x[1]]], egenes.probs[[x[2]]], 
            params)
    })
    AeqvB <- apply(pairs, 1, function(x) {
        sGenePairEqv(egenes.probs[[x[1]]], egenes.probs[[x[2]]], 
            params)
    })
    ll <- cbind(pairs, AtoB = AtoB, BtoA = BtoA, ArepB, BrepA, 
        AnonB, AeqvB)
    ll
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
