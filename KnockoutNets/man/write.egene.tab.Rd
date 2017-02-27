\name{write.egene.tab}
\alias{write.egene.tab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.egene.tab(knockdown.cols, lof, stddev, expr, file, append = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{knockdown.cols}{ ~~Describe \code{knockdown.cols} here~~ }
  \item{lof}{ ~~Describe \code{lof} here~~ }
  \item{stddev}{ ~~Describe \code{stddev} here~~ }
  \item{expr}{ ~~Describe \code{expr} here~~ }
  \item{file}{ ~~Describe \code{file} here~~ }
  \item{append}{ ~~Describe \code{append} here~~ }
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
function (knockdown.cols, lof, stddev, expr, file, append = FALSE) 
{
    if (length(knockdown.cols) != ncol(expr)) {
        stop("knockdown names don't match expression matrix dimensions")
    }
    if (length(stddev) != ncol(expr)) {
        stop("stddev length doesn't match expression matrix dimensions")
    }
    write(c("knockdown.cols:", knockdown.cols), ncolumns = length(knockdown.cols) + 
        1, sep = "\t", file = file, append = append)
    append <- TRUE
    write(c("lof:", lof), ncolumns = length(lof) + 1, sep = "\t", 
        file = file, append = append)
    write(c("stddev:", stddev), ncolumns = length(stddev) + 1, 
        sep = "\t", file = file, append = append)
    out <- cbind(rownames(expr), expr)
    colnames(out)[1] <- "Egene"
    out <- rbind(colnames(out), out)
    write.table(out, file = file, quote = FALSE, append = append, 
        sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
