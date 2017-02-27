"mixGaussians" <-
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

