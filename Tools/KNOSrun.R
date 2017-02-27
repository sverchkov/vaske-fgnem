#!/usr/bin/env Rcmd.py
library(RCommandArgs)

library(polynom)
library(logsum)
library(mixgauss)
library(lattice)
library(KnockoutNets)
library(KnockoutNetsSynth)

usage <- "

KNOSrun.R - run prediction on models generated from KNOSgen.R.  

Options:
   MODELSFILE   (required) file name of RData with models to predict
   OUTFILE      (required) file name to save predictions in R format
   TABFILE      (optional) file name to save a tab file summary
   SIGNED       (boolean) use signed data (default TRUE)
   POSMODELS    (boolean) allow positive perturbation Gaussian (default TRUE)
   TRANS        (boolean) run transitivity factors (default TRUE)
   ONLYEFFECTS  (boolean) only use true effects as data (default TRUE)
   HOLDOUT      (integer) if positive, reduce to this # SGenes (default 0)
   MARGINOP     'max' or 'sum', defaults to 'max'
"

modelsfile <- RCommandArgString("MODELSFILE", errorMsg=usage)
outfile <- RCommandArgString("OUTFILE", errorMsg=usage)
tabfile <- RCommandArgString("TABFILE", default=NA, errorMsg=usage)
signed <- RCommandArgLogical("SIGNED", default=TRUE, errorMsg=usage)
posmodels <- RCommandArgLogical("POSMODELS", default=TRUE, errorMsg=usage)
trans <- RCommandArgLogical("TRANS", default=TRUE, errorMsg=usage)
onlyeffects <- RCommandArgLogical("ONLYEFFECTS", default=T, errorMsg=usage)
holdout <- RCommandArgInteger("HOLDOUT", default=0, gte=0, errorMsg=usage)

marginop <- RCommandArgString("MARGINOP", default="max", errorMsg=usage)
marginop <- match.arg(marginop, c("max", "sum"))
marginop <- switch(marginop, max=max, sum=logsum)

loadedobjects <- load(modelsfile)

stopifnot(any("data" %in% loadedobjects))
stopifnot(class(data) == "list")

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

removeSigns <- function (x) {
  x$egenes <- abs(x$egenes);
  return(x);
}

if (!signed) {
  data <- lapply(data, removeSigns)
}

if (holdout > 0) {
  data <- lapply(data, function(d) {
    if (holdout == length(d$lof)) return(d)
    toremove <- sample(length(d$lof), length(d$lof) - holdout)
    return(removeKnockdowns(d, toremove))
  })
}

predictions <- lapply(data, function(d) {
  cat('.', file=stderr())
  params <- d$obsParams
  if (!posmodels) {
    params['pos','mean'] <- Inf
    params[,'alpha'] <- c(0.5,0.5,0)
    # must use which() here to get rid of NAs
    d$egeneParents <- d$egeneParents[which(d$egeneParents > 0)]
  }
  if (onlyeffects) {
    d$egenes <- d$egenes[which(d$egeneParents != 0), ]
  }
  return(scoreEstimates(d, params=params, doTransitivity=trans,
                        summarization=marginop))
})

save(predictions, file=outfile)

if (!is.na(tabfile)) {
  results <- collectRuns(data, predictions)
  write.table(data.frame(pair=rownames(results), results, row.names=NULL),
              file=tabfile, quote=F, sep="\t", row.names=F)  
}
