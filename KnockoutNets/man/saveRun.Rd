\name{saveRun}
\alias{saveRun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
saveRun(egeneStruct, runName, params = NULL, summarization = logsum, transiters = 100, pairscorefile = "pairscores.tab", transcorefile = "scores.tab", graphfile = "egeneobshist.png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egeneStruct}{ ~~Describe \code{egeneStruct} here~~ }
  \item{runName}{ ~~Describe \code{runName} here~~ }
  \item{params}{ ~~Describe \code{params} here~~ }
  \item{summarization}{ ~~Describe \code{summarization} here~~ }
  \item{transiters}{ ~~Describe \code{transiters} here~~ }
  \item{pairscorefile}{ ~~Describe \code{pairscorefile} here~~ }
  \item{transcorefile}{ ~~Describe \code{transcorefile} here~~ }
  \item{graphfile}{ ~~Describe \code{graphfile} here~~ }
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
function (egeneStruct, runName, params = NULL, summarization = logsum, 
    transiters = 100, pairscorefile = "pairscores.tab", transcorefile = "scores.tab", 
    graphfile = "egeneobshist.png") 
{
    if (is.null(runName) || "" == runName) {
        stop("no directory specified")
    }
    directory <- paste("Runs", runName, sep = "/")
    if (file.exists(directory)) {
        warning(paste("Directory", directory, "already exists."))
    }
    else {
        dir.create(directory, recursive = TRUE)
    }
    r <- scoreEstimates(egeneStruct, pairscorefile = pairscorefile, 
        params = params, summarization = summarization)
    appendDir <- function(f) {
        if (f == "") 
            ""
        else paste(directory, f, sep = "/")
    }
    if (graphfile != "") {
        graphfile <- appendDir(graphfile)
        png(file = graphfile, width = 800, height = 600)
        plotParameters(egeneStruct$egenes, r$params, breaks = 40)
        dev.off()
    }
    if (pairscorefile != "") {
        pairscorefile <- appendDir(pairscorefile)
        writeScoreFile(r$pairscores, pairscorefile)
    }
    if (transcorefile != "") {
        transcorefile <- appendDir(transcorefile)
        writeScoreFile(r$scores, transcorefile)
    }
    cmd <- paste("linkscores2graphviz.pl < ", sumscorefile, " > ", 
        sumscorefile, ".graphviz", sep = "")
    system(cmd)
    return(r)
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
