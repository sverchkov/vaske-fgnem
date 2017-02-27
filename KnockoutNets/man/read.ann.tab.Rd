\name{read.ann.tab}
\alias{read.ann.tab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.ann.tab(file, sep = "\t")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ ~~Describe \code{file} here~~ }
  \item{sep}{ ~~Describe \code{sep} here~~ }
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
function (file, sep = "\t") 
{
    fetchline <- function(line) {
        scan(file, what = character(0), nlines = 1, skip = line, 
            quiet = TRUE, sep = sep, quote = "")
    }
    ann <- list()
    annlines <- 0
    l <- fetchline(annlines)
    while (length(grep("^[[:alnum:].]+:$", l[1])) > 0) {
        ann[[unlist(strsplit(l[1], ":", fixed = TRUE))[1]]] <- l[2:length(l)]
        annlines <- annlines + 1
        l <- fetchline(annlines)
    }
    tab.raw <- read.delim(file = file, header = TRUE, comment.char = "", 
        quote = "", skip = annlines)
    tab <- as.matrix(tab.raw[, 2:ncol(tab.raw)])
    rownames(tab) <- tab.raw[, 1]
    c(list(tab = tab), ann)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
