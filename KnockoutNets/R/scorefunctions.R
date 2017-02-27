## ##############################
## 
## Scoring of individual edge types
##
## Breaking up each function into individual cases looks messy, but
## it's somewhat faster than the cleaner matrix-multiply
## implementation

sGenePairDomLog <- function(a.logp, b.logp, alpha,
                            summarization, na.rm=TRUE) {
  ll <- sapply(1:nrow(a.logp), function(j) {
    summarization(alpha +
                  c(a.logp[j, 'null'] + b.logp[j, 'null'], # above a
                    a.logp[j, 'neg']  + b.logp[j, 'null'], # neg, a attached
                    a.logp[j, 'neg']  + b.logp[j, 'neg'],  # neg, b attached
                    a.logp[j, 'pos']  + b.logp[j, 'null'], # pos, a attached
                    a.logp[j, 'pos']  + b.logp[j, 'pos']   # pos, b attached
                    ))
  })
  return(sum(ll, na.rm=na.rm))
}

sGenePairRepLog <- function(a.logp, b.logp, alpha,
                            summarization, na.rm=TRUE) {
  ll <- sapply(1:nrow(a.logp), function(j) {
    summarization(alpha +
                  c(a.logp[j, 'null'] + b.logp[j, 'null'], # above a
                    a.logp[j, 'neg']  + b.logp[j, 'null'], # neg, a attached
                    a.logp[j, 'pos']  + b.logp[j, 'neg'],  # neg, b attached
                    a.logp[j, 'pos']  + b.logp[j, 'null'], # pos, a attached
                    a.logp[j, 'neg']  + b.logp[j, 'pos']   # pos, b attached
                    ))
  })
  return(sum(ll, na.rm=na.rm))
}

sGenePairEqvLog <- function(a.logp, b.logp, alpha,
                            summarization, na.rm=TRUE) {
  ll <- sapply(1:nrow(a.logp), function(j) {
    summarization(alpha +
                  c(a.logp[j, 'null'] + b.logp[j, 'null'], # above a
                    a.logp[j, 'neg']  + b.logp[j, 'neg'],  # neg, b attached
                    a.logp[j, 'pos']  + b.logp[j, 'pos']   # pos, b attached
                    ))
  })
  return(sum(ll, na.rm=na.rm))
}

# negatively equivalent
sGenePairNegEqvLog <- function(a.logp, b.logp, alpha,
                            summarization, na.rm=TRUE) {
  ll <- sapply(1:nrow(a.logp), function(j) {
    summarization(alpha +
                  c(a.logp[j, 'null'] + b.logp[j, 'null'], # above a
                    a.logp[j, 'neg']  + b.logp[j, 'pos'],  # neg, b attached
                    a.logp[j, 'pos']  + b.logp[j, 'neg']   # pos, b attached
                    ))
  })
  return(sum(ll, na.rm=na.rm))
}

# Half negatively equivalent
sGenePairHalfNegEqvLog <- function(a.logp, b.logp, alpha,
                                   summarization, na.rm=TRUE) {
  ll <- sapply(1:nrow(a.logp), function(j) {
    summarization(alpha +
                  c(a.logp[j, 'null'] + b.logp[j, 'null'],
                    a.logp[j, 'neg']  + b.logp[j, 'pos'],
                    a.logp[j, 'neg']  + b.logp[j, 'neg'],
                    a.logp[j, 'pos']  + b.logp[j, 'pos'],
                    a.logp[j, 'pos']  + b.logp[j, 'neg'] 
                    ))
  })
  return(sum(ll, na.rm=na.rm))
}

sGenePairNonLog <- function(a.logp, b.logp, alpha,
                            summarization, na.rm=TRUE) {
  ll <- sapply(1:nrow(a.logp), function(j) {
    summarization(alpha +
        c(a.logp[j, 'null'] + b.logp[j, 'null'], # not attached 
          a.logp[j, 'neg']  + b.logp[j, 'null'], # neg, a attached
          a.logp[j, 'pos']  + b.logp[j, 'null'], # pos, a attached
          a.logp[j, 'null'] + b.logp[j, 'neg'],  # neg, b  attached 
          a.logp[j, 'null'] + b.logp[j, 'pos']  # pos, b  attached
        ))
  })
  return(sum(ll, na.rm=na.rm))
}

## ##############################
## 
## Prior definitions

pairPriors <- function(params, norm=thirdNullPriorNormalization) {
  return(list(dom=log(norm(domPrior(params))),
              rep=log(norm(repPrior(params))),
              eqv=log(norm(eqvPrior(params))),
              non=log(norm(nonPrior(params)))
              ))
}

# the null/null prior should be the first argument
thirdNullPriorNormalization <- function(prior, nullnull=0.8) {
  prior[1] <- 1/3
  prior[2:length(prior)] <- (2/3) / (length(prior-1)-1)
  return(prior)
  
  # prior normalized by number of regions
  return(prior / sum(prior))

  # try to get equal area under the unattached region
  prior[1] <- nullnull
  norm <- (1-nullnull) / sum(prior[2:length(prior)])
  prior[2:length(prior)] <- norm * prior[2:length(prior)]
  return(prior)
}

domPrior <- function(params) {
  prior <- c(params['null', 3] * params['null', 3],
             params['neg', 3] *  params['null', 3],
             params['neg', 3] *  params['neg', 3],
             params['pos', 3] *  params['null', 3],
             params['pos', 3] *  params['pos', 3])
  return(prior)
}

repPrior <- function(params) {
  prior <- c(params['null', 3] * params['null', 3],
             params['neg', 3]  * params['null', 3],
             params['pos', 3]  * params['neg', 3],
             params['pos', 3]  * params['null', 3],
             params['neg', 3]  * params['pos', 3])
  return(prior)
}

eqvPrior <- function(params) {
  prior <- c(params['null', 3] * params['null', 3],
             params['neg', 3] *  params['neg', 3],
             params['pos', 3] *  params['pos', 3])
  return(prior)
}

negEqvPrior <- function(params) {
  prior <- c(params['null', 3] * params['null', 3],
             params['neg', 3] *  params['pos', 3],
             params['pos', 3] *  params['neg', 3])
  return(prior)
}

halfNegEqvPrior <- function(params) {
  prior <- c(params['null', 3] * params['null', 3],
             params['neg', 3] *  params['pos', 3],
             params['neg', 3] *  params['neg', 3],
             params['pos', 3] *  params['pos', 3],
             params['pos', 3] *  params['neg', 3])
  return(prior)
}

nonPrior <- function(params) {
  prior <- c(params['null', 3] * params['null', 3],
             params['neg', 3]  * params['null', 3],
             params['pos', 3]  * params['null', 3],
             params['null', 3] * params['neg', 3],
             params['null', 3] * params['pos', 3])
  return(prior)
}

