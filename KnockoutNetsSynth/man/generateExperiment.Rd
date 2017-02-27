\name{generateExperiment}
\alias{generateExperiment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generateExperiment(nSgenes = 8, nEgenes = 20 * nSgenes, egeneInhibition = 0.25, adjGenerator = posneg, modelGenerator = randomTree, replicates = 1, params = paramGen(2, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nSgenes}{ ~~Describe \code{nSgenes} here~~ }
  \item{nEgenes}{ ~~Describe \code{nEgenes} here~~ }
  \item{egeneInhibition}{ ~~Describe \code{egeneInhibition} here~~ }
  \item{adjGenerator}{ ~~Describe \code{adjGenerator} here~~ }
  \item{modelGenerator}{ ~~Describe \code{modelGenerator} here~~ }
  \item{replicates}{ ~~Describe \code{replicates} here~~ }
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
function(nSgenes=8, nEgenes=20*nSgenes,
                               egeneInhibition=0.25,
                               adjGenerator=posneg,
                               modelGenerator=randomTree,
                               replicates=1,
                               params=paramGen(2, 1)) {

  if (is.function(nSgenes)) nSgenes <- nSgenes()
  if (is.function(nEgenes)) nEgenes <- nEgenes()
  adj <- modelGenerator(nSgenes, edgeGenerator=adjGenerator)
  acc <- knockoutEffect(adj)
  parents <- generateEgeneParents(nEgenes, colnames(adj),
                                  inhibition=egeneInhibition)
  knockdown.cols <- rep(colnames(adj), each=replicates)
  quant <- generateQuantization(knockdown.cols, parents, acc)
  expr <- generateExpression(quant, params)

  r <- list(egenes=expr, knockdown.cols=knockdown.cols, lof=colnames(adj),
            sgeneAdj=adj, sgeneAcc=acc, egeneParents=parents,
            egeneQuant=quant, obsParams=params)

  r
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
