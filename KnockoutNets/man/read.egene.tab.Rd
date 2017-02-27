\name{read.egene.tab}
\alias{read.egene.tab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.egene.tab(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (file) 
{
    result <- read.ann.tab(file)
    names(result)[1] <- "egenes"
    cols <- ncol(result$egenes)
    if (length(result$knockdown.cols) != cols) {
        stop("knockdown.cols annotation has improper length")
    }
    if (length(result$lof) == 0) {
        stop("missing lof annotation")
    }
    else {
        tmp <- rep(NA, length(result$lof))
        names(tmp) <- result$lof
        result$lof <- names(tmp)
    }
    if (!is.null(result$stddev) && length(result$stddev) != cols) {
        stop("stddev annotation length doesn't match expression matrix")
    }
    result$stddev <- as.numeric(result$stddev)
    if (any(is.na(result$stddev))) {
        stop("some stddev are not numeric")
    }
    return(result)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
