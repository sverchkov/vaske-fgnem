\name{passMessages}
\alias{passMessages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
passMessages(scoreDF, damp = 0.5, iters = 5, normToV = T, normToF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scoreDF}{ ~~Describe \code{scoreDF} here~~ }
  \item{damp}{ ~~Describe \code{damp} here~~ }
  \item{iters}{ ~~Describe \code{iters} here~~ }
  \item{normToV}{ ~~Describe \code{normToV} here~~ }
  \item{normToF}{ ~~Describe \code{normToF} here~~ }
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
function (scoreDF, damp = 0.5, iters = 5, normToV = T, normToF = NULL) 
{
    if (!isValidScoreMatrix(scoreDF)) {
        stop("invalid score matrix")
    }
    if (nrow(scoreDF) < 3) {
        return(list(scores = scoreDF))
    }
    ids <- idsOfScoreMatrix(scoreDF)
    ea <- edgeAssignments(length(ids))
    vn <- variableNeighbors(scoreDF, ea)
    constants <- as.matrix(scoreDF[, 3:ncol(scoreDF)])
    transFac <- generateTransitivityPotential()
    msgsToVars <- matrix(0, nrow = nrow(ea), ncol = 6)
    msgsToFactors <- matrix(0, nrow = nrow(ea), ncol = 6)
    oldToV <- msgsToVars
    oldToF <- msgsToFactors
    fDelta <- rep(0, iters)
    vDelta <- rep(0, iters)
    fDelta2 <- rep(0, iters)
    vDelta2 <- rep(0, iters)
    for (i in 1:iters) {
        newToF <- (1 - damp) * msgsToFactors + damp * messagesToFactor(constants, 
            msgsToVars, vn, norm = normToF)
        fDelta[i] <- sum((newToF - msgsToFactors)^2)
        fDelta2[i] <- sum((newToF - oldToF)^2)
        oldToF <- msgsToFactors
        msgsToFactors <- newToF
        newToV <- (1 - damp) * msgsToVars + damp * messagesToVars(transFac, 
            msgsToFactors, norm = normToV)
        vDelta[i] <- sum((newToV - msgsToVars)^2)
        vDelta2[i] <- sum((newToV - oldToF)^2)
        oldToV <- msgsToVars
        msgsToVars <- newToV
    }
    final <- cbind(scoreDF[, 1:2], totalVarBelief(constants, 
        msgsToVars, vn))
    colnames(final) <- colnames(scoreDF)
    result <- list(scores = final, toVars = msgsToVars, toFactors = msgsToFactors, 
        fDelta = fDelta, vDelta = vDelta, fDelta2 = fDelta2, 
        vDelta2 = vDelta2)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
