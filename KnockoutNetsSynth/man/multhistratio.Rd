\name{multhistratio}
\alias{multhistratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
multhistratio(x, breaks = "Sturges", conf.level = 0.95, length = 0.1, axislabelstart = NULL, col = gray.colors(length(x)), legend.x = NULL, legend.y = NULL, legend.txt = names(x), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{breaks}{ ~~Describe \code{breaks} here~~ }
  \item{conf.level}{ ~~Describe \code{conf.level} here~~ }
  \item{length}{ ~~Describe \code{length} here~~ }
  \item{axislabelstart}{ ~~Describe \code{axislabelstart} here~~ }
  \item{col}{ ~~Describe \code{col} here~~ }
  \item{legend.x}{ ~~Describe \code{legend.x} here~~ }
  \item{legend.y}{ ~~Describe \code{legend.y} here~~ }
  \item{legend.txt}{ ~~Describe \code{legend.txt} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function(x, breaks="Sturges",
                          conf.level=0.95, length=0.1,
                          axislabelstart=NULL,
                          col=gray.colors(length(x)),
                          legend.x=NULL, legend.y=NULL, legend.txt=names(x),
                          ...) {
  allhist <- hist(unlist(x), breaks = breaks, plot = FALSE)
  combhist <- lapply(x, function(y) {
    sapply(y, function(z) hist(z, breaks = allhist$breaks,
                                           plot = FALSE)$counts)
  })
  b <- lapply(combhist, function(x) binom.confint(x[,1], rowSums(x),
                                                  methods="exact"))
  bars <- t(sapply(b, function(x) {
    a <-  x[,'mean']
    a[x[,'n']==0] <- 0
    a
  }))
  axislabel <- signif(allhist$mids, 2)
  if (!is.null(axislabelstart)) {
    axislabel[1:length(axislabelstart)] = axislabelstart
  }
  bars.x <- barplot(bars, names=axislabel,
                    beside=TRUE, col=col, ...)
  for (i in 1:length(b)) {
    lower <- b[[i]][,'lower']
    lower[b[[i]][,'n']==0] <- NA
    upper <- b[[i]][,'upper']
    upper[b[[i]][,'n']==0] <- NA
    arrows(bars.x[i,], lower, bars.x[i,], upper,
           angle=90, code=3, length=length)
  }
  if (!is.null(legend.txt) & !is.null(legend.x) & !is.null(legend.y)) {
    legend(legend.x, legend.y, names(x), fill=col)
  }
  return(b)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
