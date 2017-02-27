\name{randomTree}
\alias{randomTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
randomTree(n, edgeGenerator = function() {
    1
}, gNames = paste("s", LETTERS[1:n], sep = ""))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{edgeGenerator}{ ~~Describe \code{edgeGenerator} here~~ }
  \item{gNames}{ ~~Describe \code{gNames} here~~ }
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
function(n, edgeGenerator=function() {1},
                       gNames=paste('s', LETTERS[1:n], sep='')) {
  r <- diag(n)
  rownames(r) <- colnames(r) <- gNames
  for (i in 2:n) {
    parent <- sample(i-1, size=1)
    r[parent, i] <- edgeGenerator()
  }
  r
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
