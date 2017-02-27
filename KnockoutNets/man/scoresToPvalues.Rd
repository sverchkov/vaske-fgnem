\name{scoresToPvalues}
\alias{scoresToPvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoresToPvalues(scores, bg.link, bg.non, bg.eqv, scorefile = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{ ~~Describe \code{scores} here~~ }
  \item{bg.link}{ ~~Describe \code{bg.link} here~~ }
  \item{bg.non}{ ~~Describe \code{bg.non} here~~ }
  \item{bg.eqv}{ ~~Describe \code{bg.eqv} here~~ }
  \item{scorefile}{ ~~Describe \code{scorefile} here~~ }
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
function (scores, bg.link, bg.non, bg.eqv, scorefile = "") 
{
    scoreSignif <- function(scores, background) {
        sapply(scores, function(x) sum(background >= x)/length(background))
    }
    s <- scores$scores
    pv <- data.frame(s[, 1:2], AtoB = scoreSignif(s[, "AtoB"], 
        bg.link), BtoA = scoreSignif(s[, "BtoA"], bg.link), AnonB = scoreSignif(s[, 
        "AnonB"], bg.non), AeqvB = scoreSignif(s[, "AeqvB"], 
        bg.eqv))
    if (scorefile != "") {
        tmp <- pv
        tmp[, 3:5] <- round(tmp[, 3:5], digits = 2)
        write.table(tmp[order(apply(s[, 3:5], 1, max), decreasing = TRUE), 
            ], file = scorefile, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    return(pv)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
