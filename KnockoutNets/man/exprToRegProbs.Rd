\name{exprToRegProbs}
\alias{exprToRegProbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
exprToRegProbs(expr, knockdown.cols, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{ ~~Describe \code{expr} here~~ }
  \item{knockdown.cols}{ ~~Describe \code{knockdown.cols} here~~ }
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
function (expr, knockdown.cols, params) 
{
    r <- lapply(unique(knockdown.cols), function(k) {
        t(apply(expr[, knockdown.cols == k, drop = FALSE], 1, 
            function(e) {
                apply(params, 1, function(d) {
                  prod(dnorm(e, mean = d[1], sd = d[2]))
                })
            }))
    })
    names(r) <- unique(knockdown.cols)
    return(r)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
