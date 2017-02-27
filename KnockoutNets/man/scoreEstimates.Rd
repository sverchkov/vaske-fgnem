\name{scoreEstimates}
\alias{scoreEstimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreEstimates(egeneStruct, summarization = logsum, pairscorefile = "", knockdown.cols = egeneStruct$knockdown.cols, egenes = egeneStruct$egenes, lof = egeneStruct$lof, runNonLof = FALSE, EMiter = 15, transIter = 100, params = NULL, doTransitivity = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egeneStruct}{ ~~Describe \code{egeneStruct} here~~ }
  \item{summarization}{ ~~Describe \code{summarization} here~~ }
  \item{pairscorefile}{ ~~Describe \code{pairscorefile} here~~ }
  \item{knockdown.cols}{ ~~Describe \code{knockdown.cols} here~~ }
  \item{egenes}{ ~~Describe \code{egenes} here~~ }
  \item{lof}{ ~~Describe \code{lof} here~~ }
  \item{runNonLof}{ ~~Describe \code{runNonLof} here~~ }
  \item{EMiter}{ ~~Describe \code{EMiter} here~~ }
  \item{transIter}{ ~~Describe \code{transIter} here~~ }
  \item{params}{ ~~Describe \code{params} here~~ }
  \item{doTransitivity}{ ~~Describe \code{doTransitivity} here~~ }
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
function (egeneStruct, summarization = logsum, pairscorefile = "", 
    knockdown.cols = egeneStruct$knockdown.cols, egenes = egeneStruct$egenes, 
    lof = egeneStruct$lof, runNonLof = FALSE, EMiter = 15, transIter = 100, 
    params = NULL, doTransitivity = T) 
{
    result <- list()
    stopifnot(typeof(logsum) == "closure")
    if (!runNonLof) {
        sub <- knockdown.cols \%in\% lof
        egenes <- egenes[, sub]
        knockdown.cols <- knockdown.cols[sub]
    }
    if (!is.null(params) && !(ncol(params) >= 3 && nrow(params) == 
        3)) {
        warning("Invalid parameter matrix, relearning")
        params <- null
    }
    if (is.null(params)) {
        params <- estimateParameters(egenes, EMiter)
    }
    result$params <- params
    kd <- unique(knockdown.cols)
    pairs <- expand.grid(A = kd, B = kd)[as.vector(upper.tri(diag(length(kd)))), 
        ]
    egenes.logprobs <- exprToRegLogProbs(egenes, knockdown.cols, 
        params)
    result$egenes.logprobs <- egenes.logprobs
    pairscores <- scorePairsWithPriorsLog(pairs, egenes.logprobs, 
        params, summarization = summarization)
    result$pairScores <- pairscores
    if (doTransitivity && (nrow(pairscores) >= 3)) {
        trans <- passMessages(pairscores, iter = transIter)
        if (trans$fDelta[transIter] != 0 || trans$vDelta[transIter] != 
            0) {
            warning("Message passing did not converge with sum scores")
        }
    }
    else {
        trans <- list(scores = pairscores)
    }
    return(c(result, trans))
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
