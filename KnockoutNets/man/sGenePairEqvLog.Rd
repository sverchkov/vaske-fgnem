\name{sGenePairEqvLog}
\alias{sGenePairEqvLog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sGenePairEqvLog(a.logp, b.logp, params, summarization, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a.logp}{ ~~Describe \code{a.logp} here~~ }
  \item{b.logp}{ ~~Describe \code{b.logp} here~~ }
  \item{params}{ ~~Describe \code{params} here~~ }
  \item{summarization}{ ~~Describe \code{summarization} here~~ }
  \item{na.rm}{ ~~Describe \code{na.rm} here~~ }
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
function (a.logp, b.logp, params, summarization, na.rm = TRUE) 
{
    prior <- log(eqvPrior(params))
    ll <- sapply(1:nrow(a.logp), function(j) {
        summarization(prior + c(a.logp[j, "null"] + b.logp[j, 
            "null"], a.logp[j, "neg"] + b.logp[j, "neg"], a.logp[j, 
            "pos"] + b.logp[j, "pos"]))
    })
    return(sum(ll, na.rm = na.rm))
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
