#!/usr/bin/env Rcmd.py

library(RCommandArgs)

library(polynom)
library(logsum)
library(mixgauss)
library(lattice)
library(KnockoutNets)

usage <- "
Calculate the likelihood of a single model, and output posteriors

    NETWORK       the network matrix
    PAIRSCORES    a pairscore file
    EGENES        tab file in egene format (defaults to stdin)
    PARAMS        the parameter matrix
    EGENESLP      an RData file with egenes log probs
    OUT           outfile
"

networkfile <- RCommandArgString("NETWORK", errorMsg=usage, default=NULL)
pairfile <- RCommandArgString("PAIRSCORES", errorMsg=usage, default=NULL)
egfile <- RCommandArgString("EGENES", default=NULL, errorMsg=usage)
paramfile <- RCommandArgString("PARAMS", default=NULL, errorMsg=usage)
eglpfile <- RCommandArgString("EGENESLP", default=NULL, errorMsg=usage)
nullnetP <- RCommandArgString("NULL", default=NULL, errorMsg=usage)
outfile <- RCommandArgFile("OUT", default="-", errorMsg=usage, open="w")

isValidParamMatrix <- function(p) {
  stopifnot(ncol(p) == 3, nrow(p) == 3)
  return(TRUE)
}

isValidNetworkMatrix <- function(n) {
  stopifnot(ncol(n) == nrow(n))
  stopifnot(!is.null(dimnames))
  stopifnot(rownames(n) == colnames(n))
  return(TRUE)
}    

if (!is.null(networkfile)) {
  network <- as.matrix(read.table(networkfile))
} else if (!is.null(pairfile)) {
  pairll <- read.delim(pairfile, stringsAsFactors=FALSE)
  stopifnot(isValidScoreMatrix(pairll))
  network <- likelihoodToAdj(pairll)
} else {
  stop("Neither PAIRSCORES nor NETWORK specified", call.=FALSE)
}
stopifnot(isValidNetworkMatrix(network))

if (is.null(eglpfile)) {
  stopifnot(!is.null(paramfile))
  params <- as.matrix(read.table(paramfile))
  stopifnot(isValidParamMatrix(params))

  stopifnot(!is.null(egfile))
  egenes <- read.egene.tab(egfile)
  # Remove non-lof genes
  sub <- egenes$knockdown.cols %in% egenes$lof
  egenes$egenes <- egenes$egenes[, sub]
  egenes$knockdown.cols <- egenes$knockdown.cols[sub]

  # "Quantize" the data
  egenes.logprobs <- exprToRegLogProbs(egenes$egenes, egenes$knockdown.cols, 
                                       params)
  stopifnot(names(egenes.logprobs) == colnames(network))
} else {
  a <- load(eglpfile)
  stopifnot("egenes.logprobs" %in% a)
}

scores <- scoreModel(network, egenes.logprobs)

cat(scores$ll, "\n")

##
##  LLR connection point
##
u <- (ncol(scores$posterior)+1)/2
connection <- colnames(scores$posterior)[-u][apply(scores$posterior[,-u],
                                               1, which.max)]
conn.sgene <- sub("^neg_", "", sub("^pos_","", connection))
conn.sign <- (regexpr("^pos_", connection) > 0) -
  (regexpr("^neg_", connection) > 0)

##
## Posterior connection strength
##
nodenames <- rownames(unique(network))
ll <- scores$posterior[,c(paste("pos_", nodenames, sep=''),
                          paste("neg_", nodenames, sep=''))]
ll <- ll + scores$egeneLL
pos <- ll - apply(ll,1,logsum)
posMax <- apply(pos, 1, max)
posRank <- rank(-posMax) # low numbers to low numbers


result <- data.frame(Egene=names(scores$maxConnectionRatio),
                     connection=conn.sgene,
                     sign=conn.sign,
                     llr=scores$maxConnectionRatio,
                     llrRank=rank(-scores$maxConnectionRatio),
                     pos=posMax,
                     posRank=posRank,
                     ll=scores$egeneLL,
                     scores$posterior)
a <- write.table(result, file=outfile, quote=FALSE, sep="\t", row.names=FALSE)
