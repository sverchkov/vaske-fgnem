\name{plotExpr}
\alias{plotExpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotExpr(egenes.probs, egeneorder = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egenes.probs}{ ~~Describe \code{egenes.probs} here~~ }
  \item{egeneorder}{ ~~Describe \code{egeneorder} here~~ }
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
function (egenes.probs, egeneorder = NULL) 
{
    lambda <- function(m) {
        val <- apply(m, 1, max)
        wh <- apply(m, 1, function(x) if (all(is.na(x))) 
            2
        else which.max(x))
        val[wh == 2] <- NA
        val[wh == 1] <- -val[wh == 1]
        val
    }
    d <- sapply(egenes.probs, lambda)
    b <- d
    b[is.na(b)] <- 0
    if (is.null(egeneorder) || length(egeneorder) != nrow(b) || 
        any(!(egeneorder \%in\% rownames(b)))) {
        egeneorder <- hclust(dist(b))$order
    }
    d <- d[egeneorder, ]
    source("~/src/cvaske/R/plotting.R")
    my.image(d, col = rgb.pal[c(1:23, 40:63)])
    return(d)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
