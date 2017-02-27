\name{plotScoreDensities}
\alias{plotScoreDensities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotScoreDensities(scoreStruct, columns, file = "", main = paste(c("Scores of ", columns)), labels = c("data", "null"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scoreStruct}{ ~~Describe \code{scoreStruct} here~~ }
  \item{columns}{ ~~Describe \code{columns} here~~ }
  \item{file}{ ~~Describe \code{file} here~~ }
  \item{main}{ ~~Describe \code{main} here~~ }
  \item{labels}{ ~~Describe \code{labels} here~~ }
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
function (scoreStruct, columns, file = "", main = paste(c("Scores of ", 
    columns)), labels = c("data", "null")) 
{
    null <- scoreStruct$scores[scoreStruct$lofVSnonlof, columns]
    data <- scoreStruct$scores[scoreStruct$lofVSlof, columns]
    null.dens <- density(unlist(null))
    data.dens <- density(unlist(data))
    if (file != "") {
        png(file = file)
    }
    plotDensities(data.dens, null.dens, main = main, labels = labels)
    if (file != "") {
        dev.off()
    }
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
