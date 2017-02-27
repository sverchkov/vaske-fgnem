\name{norm.kl.div}
\alias{norm.kl.div}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
norm.kl.div(m1, s1, m2, s2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{ ~~Describe \code{m1} here~~ }
  \item{s1}{ ~~Describe \code{s1} here~~ }
  \item{m2}{ ~~Describe \code{m2} here~~ }
  \item{s2}{ ~~Describe \code{s2} here~~ }
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
function (m1, s1, m2, s2) 
{
    return(0.5 * (log(s2/s1) + s1/s2 + (m1 - m2) * (m1 - m2)/s2 - 
        1))
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
