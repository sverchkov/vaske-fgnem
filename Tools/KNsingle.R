#!/usr/bin/env Rcmd.py
library(RCommandArgs)
library(Cairo)

library(polynom)
library(logsum)
library(mixgauss)
library(lattice)
library(KnockoutNets)

main <- function() {
  usage <- "
Score a single network, printing the scores to std out.

Options:
    EGENES        tab file of Egenes
    MEAN          offset of the mean for the up and down observation
                      distribution (default is 1.5)
    SD            Standard deviation of observation distributions (default 1)
    ESTIMATEP     Run EM to estimate parameters. If set, either 'balanced'
                  or 'unbalanced'
    POSONLY       restrict to only positive connections (default FALSE)
    EGENESLP      File of estimated egenes.logprobs
    MARGINOP     'max' or 'sum', defaults to 'max'
    DOPART        (optional) which part of the subset of pairs to score
    NUMPARTS      (optional) the number of parts to score
    TRANSITIVE    do transitive message passing (default FALSE)
    SAVERESULTS   if set, save an R object of the results to this file
    PERMUTATIONS  number of permutations to perform
    NETWORKSCORE  where to save the results of entire network scoring
    PLOTOBSHIST   if set, make a histogram of the observations with the
                  parameters.  The value should be an R expression
                  that opens a device.
    PLOTOBSIMAGE  if set, make a heat map of the observations with the
                  likelihood, colored by the distribution from which that
                  likelihood came.  The value shoud be an R expression
"

  egeneFile <- RCommandArgString("EGENES", default=NULL, errorMsg=usage)
  paramMean <- RCommandArgDouble("MEAN", default=1.5, gt=0, errorMsg=usage)
  paramSD <- RCommandArgDouble("SD", default=1, gt=0)
  estimateP <- RCommandArgString("ESTIMATEP", default=NA)
  posonly <- RCommandArgLogical("POSONLY", default=FALSE)

  egenes.lpfile <- RCommandArgString("EGENESLP", default=NULL, errorMsg=usage)

  stopifnot(xor(is.null(egeneFile), is.null(egenes.lpfile)))

  marginop <- RCommandArgString("MARGINOP", default="max", errorMsg=usage)
  marginop <- match.arg(marginop, c("max", "sum"))
  marginop <- switch(marginop, max=max, sum=logsum)

# Haven't started on this yet
                                        #numParts <- RCommandArgInteger("NUMPARTS", gte=1, default=1 errorMsg=usage)
                                        #doPart <- RCommandArgInteger("NUMPARTS", gte=1, lte=numParts,
                                        #                             default=1 errorMsg=usage)

  trans <- RCommandArgLogical("TRANSITIVE", default=FALSE)
  rsave <- RCommandArgString("SAVERESULTS", default=NULL)
  perms <- RCommandArgInteger("PERMUTATIONS", gte=0, default=0)
  netscorefile <- RCommandArgString("NETWORKSCORE", default=NULL)
  plotobshist <- RCommandArgString("PLOTOBSHIST", default=NULL)
  plotobsimage <- RCommandArgString("PLOTOBSIMAGE", default=NULL)

  if(!is.null(egeneFile)) {
    eg <- read.egene.tab(egeneFile)
    params <- paramGen(paramMean, paramSD)
    if (!is.na(estimateP)) {
      estimateP <- pmatch(estimateP, c("unbalanced", "balanced"))[1]
      if (is.na(estimateP)) {
        stop("unknown parameter estimating method")
      }
      x <- eg$egenes[,eg$knockdown.cols %in% eg$lof]
      if (estimateP == 2) {
        x <- c(x,-x)
      }
      params <- estimateParameters(x, mu=params[,'mean'], sigma=params[,'sd'])
    }
    
    if (posonly) {
      params['pos','mean'] <- Inf
      posalpha <- params['pos','alpha']
      params[,'alpha'] <- params[,'alpha'] +
        c(posalpha/2, posalpha/2, -posalpha)
    }
  } else {
    stop("EGenes.LP part not implemented")
  }
    
                                        # Haven't started on this yet
                                        # doParts <- NULL
                                        # if (numParts > 1) {
                                        #   doParts <- c(doPart, numParts)
                                        # }
  results <- scoreBestModelEstimate(eg, params=params, doTransitivity=trans,
                                        #                                  doPart=doParts, 
                                    summarization=marginop)

  a <- write.table(results$scores, file=stdout(), quote=F, row.names=F, sep="\t")

  permRun <- function() {
    e <- eg;
    lof <- e$knockdown.cols %in% e$lof
    #e$egenes[, lof] <- sample(e$egenes[,lof])
    e$egenes[,lof] <- apply(e$egenes[,lof], 1, sample)
    return(scoreBestModelEstimate(e, params=params, doTransitivity=trans,
                                  summarization=marginop))
  }
  permutations <- replicate(perms, permRun(), simplify=FALSE)

  if(!is.null(rsave)) {
    save(results, file=rsave)
  }

  if(!is.null(netscorefile) && perms > 0) {
    write(c(results$ll, sapply(permutations, function(x) x$ll)), sep="\t",
          ncolumns=(length(permutations)+1), file=netscorefile)
  }

  if(!is.null(plotobshist)) {
    eval(parse(text=plotobshist))
    plotParameters(eg$egenes[,eg$knockdown.cols %in% eg$lof], params)
    a <- dev.off()
  }
  if(!is.null(plotobsimage)) {
    eval(parse(text=plotobsimage))
    plotExpr(results$egenes.logprobs)
    a <- dev.off()
  }
}

#debug(main)

main()
