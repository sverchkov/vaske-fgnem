\name{mixGaussians}
\alias{mixGaussians}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ mixGaussians }
\description{
Basic mixture of Gaussians algorithm.  Given a set of one 
  dimensional points and a number of Gaussians, use EM to determine each 
  Gaussian's mean and standard deviation.
}
\usage{
mixGaussians(x, models, iter.max = 15, mu = sample(unique(x), models), sigma = rep(sd(x)/models, models), fixed.mu = rep(NA, models), fixed.alpha = rep(NA, models))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{models}{ ~~Describe \code{models} here~~ }
  \item{iter.max}{ ~~Describe \code{iter.max} here~~ }
  \item{mu}{ ~~Describe \code{mu} here~~ }
  \item{sigma}{ ~~Describe \code{sigma} here~~ }
  \item{fixed.mu}{ ~~Describe \code{fixed.mu} here~~ }
  \item{fixed.alpha}{ ~~Describe \code{fixed.alpha} here~~ }
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
function(x, models, iter.max=15,
                         mu=sample(unique(x), models),
                         sigma=rep(sd(x)/models, models),
                         fixed.mu=rep(NA, models),
                         fixed.alpha=rep(NA, models)) {
  if (length(models) != 1 || models != as.integer(models)) {
    stop("models should be a single integer")
  }

  # mixture parameters
  alpha <- rep(1/models, models)

  # mean for each model
  if (is.null(mu) || length(mu) != models) {
    mu <- sample(unique(x), models)
  }
  mu[!is.na(fixed.mu)] <- fixed.mu[!is.na(fixed.mu)]

  # standard deviation for each model
  if (is.null(sigma) || length(sigma) != models) {
    sigma <- rep(sd(x)/models, models)
  }

                                        # each row is P(y_i | x_i, THETA)
  y <- matrix(1/models, nrow=length(x), ncol=models)

  for (iter in 1:iter.max) {
    ##
    ## Estimation
    ##
    ynew <- sapply(1:models,
                   function(i) alpha[i] * dnorm(x, mean=mu[i], sd=sigma[i]))
    weight <- apply(ynew, 1, sum)
    ynew <- ynew / weight
    y <- ynew
    modelWeights <- colSums(y)

    ##
    ## Maximization
    ##
    
    alpha <- colMeans(y)

    mu <- colSums(y * x) / modelWeights
    mu[!is.na(fixed.mu)] <- fixed.mu[!is.na(fixed.mu)]

    sigma <- sqrt(sapply(1:models, function(m) {
      d <- x-mu[m]
      sum(d*d*y[,m]) / modelWeights[m]
    }))
  }

  result <- list(y=y, mu=mu, sigma=sigma, alpha=alpha)
  class(result) <- "mixGaussians"
  result
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
