\name{scoreModel}
\alias{scoreModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreModel(accMatrix, egenes.logprobs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{accMatrix}{ ~~Describe \code{accMatrix} here~~ }
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
function (accMatrix, egenes.logprobs) 
{
    stopifnot(rownames(accMatrix) == names(egenes.logprobs), 
        colnames(accMatrix) == names(egenes.logprobs))
    logprobs <- log.probsTo3Array(egenes.logprobs)
    diag(accMatrix) <- 1
    upMatrix <- accMatrix + 2
    downMatrix <- 2 - accMatrix
    pickFromRow <- function(m, pick) {
        m[0:(ncol(m) - 1) * nrow(m) + pick]
    }
    getParentLikelihoods <- function(distByKO, accToDistMapping) {
        apply(accToDistMapping, 2, function(KOEffect) {
            sum(pickFromRow(distByKO, KOEffect))
        })
    }
    al <- apply(logprobs, 1, function(distByKO) {
        c(getParentLikelihoods(distByKO, upMatrix), sum(distByKO[2, 
            ]), getParentLikelihoods(distByKO, downMatrix))
    })
    rownames(al) <- c(paste("neg", rownames(accMatrix), sep = "_"), 
        "unattached", paste("pos", rownames(accMatrix), sep = "_"))
    zs <- apply(al, 2, function(x) logsum(sort(x)))
    posterior <- t(al) - zs
    ll <- sum(zs, na.rm = T)
    return(list(ll = ll, posterior = posterior))
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
