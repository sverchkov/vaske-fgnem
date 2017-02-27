\name{writeScoreFile}
\alias{writeScoreFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeScoreFile(scores, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{ ~~Describe \code{scores} here~~ }
  \item{file}{ ~~Describe \code{file} here~~ }
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
function (scores, file) 
{
    tmp <- scores
    tmp[, 3:8] <- round(tmp[, 3:8], digits = 4)
    write.table(tmp[order(apply(scores[, 3:8], 1, max), decreasing = TRUE), 
        ], file = file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
