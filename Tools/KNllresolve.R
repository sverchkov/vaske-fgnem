#!/usr/bin/env Rgetopt.py

library(Rgetopt)
library(logsum)
library(mixgauss)
library(KnockoutNets)

argspec <- c("\
KNllresolve.R - read in a pairwise score matrix, do transitivity resolution,
perform transitive reduction ala Wagner 2001, and output result.",
             "data|d=rfile     Pairwise score matrix",
             "iters|i=i        iterations of message passing (default 50)",
             "llfile=wfile     If set, write ll score file here",
             "noreduce         Skip the transitive reduction",
             "dotfile=wfile    If set, write a dot file here")
o <- Rgetopt(argspec)

##
## Check options
##

if (is.null(o$data)) usage("Must specify a data file", argspec)
if (length(o$argv) > 0) usage(paste("extra unparsed arguments:",
                                    o$argv, collapse=" "), argspec)
iters <- if(is.null(o$iters)) 50 else o$iters
stopifnot(iters >= 0)

if (is.null(o$dotfile) && is.null(o$llfile)) {
  stop("No output files specified")
}
             
##
## read, process
##

ll <- read.delim(o$data, header=T, row.names=NULL)

stopifnot(isValidScoreMatrix(ll))

trans <- passMessages(ll, iters=iters)

if (!is.null(o$llfile)) {
  write.table(trans$scores, file=o$llfile, quote=F, row.names=F, sep="\t")
}

if (!is.null(o$dotfile)) {
  net <- likelihoodToAdj(trans$scores)
  
  reduceError <- NULL
  net <- mergeEquiv(net)
  if (is.null(o$noreduce)) {
    net <- mergeEquiv(signedAcc(net))
    reduction <- try(wagnerReduction(net))
    if ("try-error" %in% class(reduction)) {
      reduceError <- c("#transitive reduction failed, perhaps a cycle?")
    } else {
      net <- reduction
    }
  }

  dot <- accToDot(net)
  writeLines(c(reduceError, dot), con=o$dotfile)
}
