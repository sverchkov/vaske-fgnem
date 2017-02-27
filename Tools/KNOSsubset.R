#!/usr/bin/env Rgetopt.py
library(Rgetopt)

argspec <- c("KNOSsubset.R - ",
             "modelsfile=s   file name of RData with models to predict",
             "holdout=i      number of Sgenes to keep",
             "o=s            where to write output file")

removeKnockdowns <- function(d, toremove) {
  hybremove <- which(d$knockdown.cols %in% d$lof[toremove])
  d$egenes <- d$egenes[,-hybremove]
  d$knockdown.cols <- d$knockdown.cols[-hybremove]
  d$lof <- d$lof[-toremove]
  d$sgeneAdj <- d$sgeneAdj[-toremove, -toremove]
  d$sgeneAcc <- d$sgeneAcc[-toremove, -toremove]
  d$egeneParents[d$egeneParents %in% toremove] <- NA
  d$egeneQuant <- d$egeneQuant[,-hybremove]
  # obsParams is fine
  return(d)
}

main <- function(argv) {
  if (missing(argv)) argv <- RgetArgvEnvironment()[-1]

  o <- Rgetopt(argv=argv, argspec=argspec)

  stopifnot(!is.null(o$o), !is.null(o$modelsfile), !is.null(o$holdout))

  loadedobjects <- load(o$modelsfile)

  if (o$holdout > 0) {
    data <- lapply(data, function(d) {
      if (o$holdout == length(d$lof)) return(d)
      toremove <- sample(length(d$lof), length(d$lof) - o$holdout)
      return(removeKnockdowns(d, toremove))
    })
  }

  save(data, file=o$o)
}

main()
