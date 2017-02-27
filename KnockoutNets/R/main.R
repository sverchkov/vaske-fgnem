#' Write the E-gene tab-separated file
#'
#' @param knockdown.cols Vector of knockdown names (one per expression matrix column)
#' @param lof Loss-of-function genes (vector)
#' @param stddev Vector of standard deviations (one per expression matrix column)
#' @param expr Matrix of gene expression levels
#' @param file Output file
#' @param append Whether to append to the file
#' @export
write.egene.tab <- function(knockdown.cols, lof, stddev, expr,
                            file, append=FALSE) {
  if (length(knockdown.cols) != ncol(expr)) {
    stop ("knockdown names don't match expression matrix dimensions")
  }
  if (length(stddev) != ncol(expr)) {
    stop ("stddev length doesn't match expression matrix dimensions")
  }

  write(c("knockdown.cols:", knockdown.cols),
        ncolumns=length(knockdown.cols)+1,
        sep="\t", file=file, append=append)
  append <- TRUE

  write(c("lof:", lof), ncolumns=length(lof)+1,
        sep="\t", file=file, append=append)

  write(c("stddev:", stddev), ncolumns=length(stddev)+1,
        sep="\t", file=file, append=append)

  out <- cbind(rownames(expr), expr)
  colnames(out)[1] <- "Egene"
  out <- rbind(colnames(out), out)
  write.table(out, file=file, quote=FALSE, append=append, sep="\t",
              row.names=FALSE, col.names=FALSE)
}

#' @export
read.ann.tab <- function(file, sep="\t") {
  fetchline <- function(line) {
    scan(file, what=character(0), nlines=1, skip=line, quiet=TRUE, sep=sep,
         quote="")
  }

  ann <- list()
  annlines <- 0
  l <- fetchline(annlines)
  while(length(grep('^[[:alnum:].]+:$', l[1])) > 0) {
    ann[[unlist(strsplit(l[1], ':', fixed=TRUE))[1]]] <- l[2:length(l)]
    annlines <- annlines + 1
    l <- fetchline(annlines)
  }
  tab.raw <- read.delim(file=file, header=TRUE, comment.char="",
                    quote="", skip=annlines)
  tab <- as.matrix(tab.raw[,2:ncol(tab.raw)])
  rownames(tab) <- tab.raw[,1]
  
  c(list(tab=tab), ann)
}

#' @export
read.egene.tab <- function(file) {
  result <- read.ann.tab(file)
  names(result)[1] <- "egenes"
  cols <- ncol(result$egenes)
  if (length(result$knockdown.cols) != cols) {
    stop("knockdown.cols annotation has improper length")
  }
  if (length(result$lof) == 0) {
    stop("missing lof annotation")
  } else {
    tmp <- rep(NA, length(result$lof))
    names(tmp) <- result$lof
    result$lof <- names(tmp)
  }
  if (!is.null(result$stddev) && length(result$stddev) != cols) {
    stop("stddev annotation length doesn't match expression matrix")
  }
  result$stddev <- as.numeric(result$stddev)
  if (any(is.na(result$stddev))) {
    stop("some stddev are not numeric")
  }
  return(result)
}

#' @export
paramGen <- function(offset, sd) {
  matrix(c(-offset, 0, offset, sd, sd, sd, 1/3, 1/3, 1/3), nrow=3, ncol=3,
         dimnames=list(c("neg", "null", "pos"), c("mean", "sd", "alpha")))
}

#' Convert from expression data to scores
#' params should be a vector of sds of length equal to number of columns in expr
#' @export
exprToRegLogProbs <- function(expr, knockdown.cols, params, summarize=sum) {
  r <- lapply(unique(knockdown.cols), function(k) {
    t(apply(expr[,knockdown.cols==k, drop=FALSE], 1, function(e) {
      apply(params, 1, function(d) {
        summarize(dnorm(e, mean=d[1], sd=d[2], log=TRUE)) +
          length(e) * log(d[3])
      })
    }))
  })
  names(r) <- unique(knockdown.cols)
  return(r)
}

#' Convert from expression data to log scores
#' params should be a matrix 3x(2+), with the first column being the means
#' and the second column being the standard deviations
#' @export
exprToRegProbs <- function(expr, knockdown.cols, params) {
  n <- unique(knockdown.cols)
  r <- lapply(n, function(k) {
    t(apply(expr[,knockdown.cols==k, drop=FALSE], 1, function(e) {
      apply(params, 1, function(d) {
        prod(dnorm(e, mean=d[1], sd=d[2]))
      })
    }))
  })
  names(r) <- n
  return(r)
}

#' @export
estimateParameters <- function(egenes, EMiter=15, mu=NULL, sigma=NULL,
                               verbose=FALSE) {
  x <- as.vector(egenes)
  x <- x[!is.na(x)]
  mg <- mixGaussians(x, models=3, iter.max=EMiter,
                     mu=mu, sigma=sigma, fixed.mu=c(NA, 0, NA))

  if (verbose) {
    print(mg)
  }
  
  params <- cbind(mg$mu, mg$sigma, mg$alpha)
  
  params <- params[order(mg$mu),]
  colnames(params) <- c("mean", "sd", "alpha")
  rownames(params) <- c("neg", "null", "pos")
  return(params)
}

#' @export
estimateParameters2 <- function(egenes, EMiter=15, mu=NULL, sigma=NULL) {
  x <- as.vector(egenes)
  x <- x[!is.na(x)]
  mg <- mixGaussians(x, models=2, iter.max=EMiter,
                     mu=mu, sigma=sigma, fixed.mu=c(0, NA))
    
  params <- cbind(mg$mu, mg$sigma, mg$alpha)
  
  colnames(params) <- c("mean", "sd", "alpha")
  if (sum(params[,'mean']) > 0) {
    params <- rbind(c(Inf, 0, 0), params)
  } else {
    params <- rbind(c(Inf, 0, 0), params)[3:1,]
  }
  rownames(params) <- c("neg", "null", "pos")
  return(params)
}

##' Scoring using true probabilities.  I've found this to be less accurate
##' then using the log probabilities (scorePairsWithPriorsLog)
##'
##' @export
scorePairsWithPriors <- function(pairs, egenes.probs, egenes.logprobs,
                                 params) {
  AtoB <- apply(pairs, 1, function(x) {
    sGenePairDom(egenes.probs[[x[1]]], egenes.probs[[x[2]]], params)
  })
  BtoA <- apply(pairs, 1, function(x) {
    sGenePairDom(egenes.probs[[x[2]]], egenes.probs[[x[1]]], params)
  })
  ArepB <- apply(pairs, 1, function(x) {
    sGenePairRep(egenes.probs[[x[1]]], egenes.probs[[x[2]]], params)
  })
  BrepA <- apply(pairs, 1, function(x) {
    sGenePairRep(egenes.probs[[x[2]]], egenes.probs[[x[1]]], params)
  })
  AnonB <- apply(pairs, 1, function(x) {
    sGenePairNon(egenes.probs[[x[1]]], egenes.probs[[x[2]]], params)
  })
  AeqvB <- apply(pairs, 1, function(x) {
    sGenePairEqv(egenes.probs[[x[1]]], egenes.probs[[x[2]]], params)
  })
S
  ll <- cbind(pairs, AtoB=AtoB, BtoA=BtoA, ArepB, BrepA, AnonB, AeqvB)
  ll
}

#' @export
scorePairsWithPriorsLog <- function(pairs, egenes.logprobs,
                                    alphas, summarization) {
  AtoB <- apply(pairs, 1, function(x) {
    sGenePairDomLog(egenes.logprobs[[x[1]]], egenes.logprobs[[x[2]]],
                    alpha=alphas$dom, summarization=summarization)
  })
  BtoA <- apply(pairs, 1, function(x) {
    sGenePairDomLog(egenes.logprobs[[x[2]]], egenes.logprobs[[x[1]]],
                    alpha=alphas$dom, summarization=summarization)
  })
  ArepB <- apply(pairs, 1, function(x) {
    sGenePairRepLog(egenes.logprobs[[x[1]]], egenes.logprobs[[x[2]]],
                    alpha=alphas$rep, summarization=summarization)
  })
  BrepA <- apply(pairs, 1, function(x) {
    sGenePairRepLog(egenes.logprobs[[x[2]]], egenes.logprobs[[x[1]]],
                    alpha=alphas$rep, summarization=summarization)
  })
  AnonB <- apply(pairs, 1, function(x) {
    sGenePairNonLog(egenes.logprobs[[x[1]]], egenes.logprobs[[x[2]]],
                    alpha=alphas$non, summarization=summarization)
  })
  AeqvB <- apply(pairs, 1, function(x) {
    sGenePairEqvLog(egenes.logprobs[[x[1]]], egenes.logprobs[[x[2]]],
                    alpha=alphas$eqv, summarization=summarization)
  })

  ll <- cbind(pairs, AtoB=AtoB, BtoA=BtoA, ArepB, BrepA, AnonB, AeqvB)
  ll
}

#' @export
saveRun <- function(egeneStruct, runName,
                    params=NULL,
                    summarization=logsum,
                    transiters=100,
                    pairscorefile="pairscores.tab",
                    transcorefile="scores.tab",
                    graphfile="egeneobshist.png") {
  ## Prepare the output directory
  if (is.null(runName) || "" == runName) {
    stop("no directory specified")
  }
  directory <- paste("Runs", runName, sep="/")
  if (file.exists(directory)) {
    warning(paste("Directory", directory, "already exists."))
  } else {
    dir.create(directory, recursive=TRUE)
  }
  
  ## Do the run
  r <- scoreEstimates(egeneStruct,
                      pairscorefile=pairscorefile,
                      params=params,
                      summarization=summarization)

  # Output all the data
  appendDir <- function(f) {if (f == '') '' else paste(directory, f, sep='/')}
  if (graphfile != '') {
    graphfile <- appendDir(graphfile)
    png(file=graphfile, width=800, height=600)
    plotParameters(egeneStruct$egenes, r$params, breaks=40)
    dev.off()
  }
  if (pairscorefile != '') {
    pairscorefile <- appendDir(pairscorefile)
    writeScoreFile(r$pairscores, pairscorefile)
  }
  if (transcorefile != '') {
    transcorefile <- appendDir(transcorefile)
    writeScoreFile(r$scores, transcorefile)
  }

  cmd <- paste("linkscores2graphviz.pl < ", sumscorefile,
               " > ", sumscorefile, ".graphviz", sep='')
  system(cmd)

  return(r)
}

#' @export
scoreEstimates <- function(egeneStruct, 
                           summarization=logsum,
                           pairscorefile='',
                           knockdown.cols=egeneStruct$knockdown.cols,
                           egenes=egeneStruct$egenes,
                           lof=egeneStruct$lof,
                           runNonLof=FALSE,
                           EMiter=15,
                           transIter=100,
                           params=NULL,
                           alphaPriors=pairPriors(params),
                           doPart=NULL,
                           doTransitivity=T) {
  result <- list()

  summarization <- match.fun(summarization)

  # Subset if we're not using arrays that exhibit the phenotype
  if (!runNonLof) {
    sub <- knockdown.cols %in% lof
    egenes <- egenes[,sub]
    knockdown.cols <- knockdown.cols[sub]
  }

  # Ensure we have parameters
  if (!is.null(params) && !(ncol(params) >= 3 && nrow(params) == 3)) {
    warning("Invalid parameter matrix, relearning")
    params <- null
  }
  if (is.null(params)) {
    params <- estimateParameters(egenes, EMiter)
  }
  result$params <- params

  # Get all the pairs
  kd <- unique(knockdown.cols)
  pairs <- expand.grid(A=kd,B=kd)[as.vector(upper.tri(diag(length(kd)))),]

  if (!is.null(doPart)) {
    stopifnot(length(doPart) == 2)
    stopifnot(doPart[1] >= 1, doPart[1] <= doPart[2])

    f <- rep(1:doPart[2], length.out=nrow(pairs))
    pairs <- pairs[f == doPart[1], ]
    kd <- unique(unlist(pairs))
    s <- knockdown.cols %in% kd
    egenes <- egenes[,s]
    knockdown.cols <- knockdown.cols[,s]
  }
  
  # Quantize the data
  egenes.logprobs <- exprToRegLogProbs(egenes, knockdown.cols, params)
  result$egenes.logprobs <- egenes.logprobs

  # Estimate the pairwise scores
  pairscores <- scorePairsWithPriorsLog(pairs, egenes.logprobs,
                                        alphas=alphaPriors,
                                        summarization=summarization)
  result$pairScores <- pairscores

  # Resolve transitivity
  if (doTransitivity && (nrow(pairscores) >= 3) && is.null(doPart)) {
    trans <- passMessages(pairscores, iter=transIter)
    if (trans$fDelta[transIter] != 0 || trans$vDelta[transIter] != 0) {
      warning("Message passing did not converge with sum scores")
    }
  } else {
    trans <- list(scores=pairscores)
  }

  return(c(result, trans))
}
  
#' take an n x 2 character matrix, each row representing a link
#' and convert to an adjacency matrix
#' @export
generateAdj <- function(links) {
  Sgenes <- unique(as.vector(links))
  adj <- diag(length(Sgenes))
  rownames(adj) <- colnames(adj) <- Sgenes
  apply(links, 1, function(x) adj[x[1], x[2]] <<- 1)
  return(adj)
}

#' Thisk looks like a very inefficient way to do this...
#' @export
matPower <- function(X,n){
    if(n != round(n)) {
        n <- round(n)
        warning("rounding exponent `n' to", n)
    }
    phi <- diag(nrow = nrow(X))
    pot <- X # the first power of the matrix.

    while (n > 0)
    {
        if (n %% 2)
            phi <- phi %*% pot

        n <- n %/% 2
        pot <- pot %*% pot
    }
    return(phi)
}

#' @export
scoresToPvalues <- function(scores, bg.link, bg.non, bg.eqv, scorefile='') {
  scoreSignif <- function(scores, background) {
    sapply(scores, function(x)  sum(background >= x)/length(background))
  }
  s <- scores$scores
  pv <- data.frame(s[,1:2],
                   AtoB=scoreSignif(s[,"AtoB"], bg.link),
                   BtoA=scoreSignif(s[,"BtoA"], bg.link),
                   AnonB=scoreSignif(s[,"AnonB"], bg.non),
                   AeqvB=scoreSignif(s[,"AeqvB"], bg.eqv))
  if (scorefile != '') {
    tmp <- pv
    tmp[,3:5] <- round(tmp[,3:5], digits=2)
    # order by score, not pvalue (as pvalue loses ordering when tied)
    write.table(tmp[order(apply(s[,3:5],1,max), decreasing=TRUE),], 
                file=scorefile, sep="\t", quote=FALSE, row.names=FALSE)
  }
  return(pv)
}

## ##############################
##
## Output functions
##

#' @export
plotParameters <- function(egenes, p, ...) {
  egenes <- egenes[!is.na(egenes)]
  mesh <- seq(from=range(egenes)[1], to=range(egenes)[2], length.out=500)
  hist(as.vector(egenes), prob=TRUE, main='Egene observations', ...)
  linecol <- c(3,1,2)[rank(p[,'mean'])]
  for (i in 1:nrow(p)) {
    lines(mesh, dnorm(mesh, mean=p[i,'mean'], sd=p[i,'sd']) * p[i,'alpha'],
          col=linecol[i])
  }
}

#' @export
writeScoreFile <- function(scores, file) {
    tmp <- scores
    tmp[,3:8] <- round(tmp[,3:8], digits=4)
    write.table(tmp[order(apply(scores[,3:8],1,max), decreasing=TRUE),],
                file=file, sep="\t", quote=FALSE, row.names=FALSE)
}

#' @export
plotScoreDensities <- function(scoreStruct, columns, file='',
                               main=paste(c('Scores of ', columns)),
                               labels=c("data", "null")) {
  null <- scoreStruct$scores[scoreStruct$lofVSnonlof, columns]
  data <- scoreStruct$scores[scoreStruct$lofVSlof, columns]
  
  null.dens <- density(unlist(null))
  data.dens <- density(unlist(data))

  if (file != '') {
    png(file=file)
  }

  plotDensities(data.dens, null.dens, main=main, labels=labels)
  
  if (file != '') {
    dev.off()
  }
}

#' @export
plotDensities <- function (densa, densb, acol=1, bcol=2,
                          labels=character(0), ...) {
  xlim = range(c(densa$x, densb$x))
  ylim = range(c(densa$y, densb$y))

  plot(densa, xlim=xlim, ylim=ylim, col=acol, ...)
  lines(densb$x, densb$y, col=bcol)
  if (length(labels) == 2) {
    legend(x=xlim[1], y=ylim[2], legend=labels, col=c(acol,bcol), lty=1)
  }
}

#' @export
revEgeneStruct <- function(egeneStruct) {
  r <- egeneStruct
  r$stddev <- rev(r$stddev)
  r$knockdown.cols <- rev(r$knockdown.cols)
  r$egenes <- r$egenes[,ncol(r$egenes):1]
  r
}

#' @export
plotProbs <- function(scoreStruct, scaleLiklihood=F, zlim=c(-5,5)) {
  m <- sapply(scoreStruct$egenes.probs, function(x) x[,3]-x[,1])
  if (scaleLiklihood) {
    r <- range(m)
    m[m<0] <- m[m<0] / -r[1]
    m[m>0] <- m[m>0] / r[2]
  }
  if (!is.null(zlim)) {
    m[m<zlim[1]] <- zlim[1]
    m[m>zlim[2]] <- zlim[2]
  }
  source("~/src/cvaske/R/plotting.R")
  yb.pal <- c(rgb(0,0,16:0/16), rgb(1:16/16, 1:16/16,0))
  my.image(m, zlim=zlim, col=yb.pal)
}

#' @export
plotExpr <- function(egenes.probs, egeneorder=NULL, doplot=T) {
  lambda <- function (m) {
    val <- apply(m, 1, max)
    wh <- apply(m, 1, function (x) if (all(is.na(x))) 2 else which.max(x))
    val[wh==2] <- NA
    val[wh==1] <- -val[wh==1]
    val
  }
  d <-sapply(egenes.probs, lambda)
  b <- d
  b[is.na(b)] <- 0
  if (is.null(egeneorder) || length(egeneorder) != nrow(b)
      || any(!(egeneorder %in% rownames(b)))) {
    egeneorder <- hclust(dist(b))$order
  }
  d <- d[egeneorder,]
  if(doplot) {
    source("~/src/cvaske/R/plotting.R")
    my.image(d, col=rgb.pal[c(1:23,40:63)])
  }
  return(d)
}

#' @export
norm.kl.div <- function(m1, s1, m2, s2) {
  return(.5 * (log(s2/s1) + s1/s2 + (m1-m2)*(m1-m2)/s2 - 1))
}

#' @export
param.div <- function(params) {
  outer(rownames(params), rownames(params), FUN=function(x,y) {
    norm.kl.div(params[x,'mean'], params[x,'sd'],
                params[y,'mean'], params[x,'sd'])
  })
}

#' @export
signedAcc <- function(adjMatrix) {
  if (!(ncol(adjMatrix) > 0) || ncol(adjMatrix) != nrow(adjMatrix)) {
    stop("Bad input matrix")
  }

  diag(adjMatrix) <- 0
  acc <- adjMatrix %*% diag(1, nrow=nrow(adjMatrix))
  for (i in 1:nrow(adjMatrix)) {
    acc0 <- acc == 0
    accNew <- (acc %*% adjMatrix)[acc0]
    if (any(accNew != 0)) {
      acc[acc0] <- sign(accNew)
    } else {
      break
    }
  }
  diag(acc) <- 1
  colnames(acc) <- rownames(acc)
  return(acc)
}

#' Presume that that everything upstream/out of stream of the
#' knockdown is activated, and follow through with knockdown
#' expression of the
#' @export
knockoutEffect <- function(adjMatrix) {
  acc <- matrix(0, nrow=nrow(adjMatrix), ncol=nrow(adjMatrix),
                dimnames=dimnames(adjMatrix))
  accNew <- adjMatrix
  while (any((acc == 0) & (accNew != 0))) {
    better <- abs(accNew) > abs(acc)
    acc[better] <- accNew[better]
    accNew <- MatMultAbsMax(acc, adjMatrix)
  }
  return(acc)
}

#' @export
MatMultAbsMax <- function(A, B) {
  t(apply(A, 1, function(x) {
    apply(x*B, 2, function(y) y[which.max(abs(y))])
  }))
}

#' @export
log.probsTo3Array <- function(egenes.logprobs) {
  array(unlist(egenes.logprobs),
        c(dim(egenes.logprobs[[1]]), length(egenes.logprobs)),
        dimnames=list(rownames(egenes.logprobs[[1]]),
          colnames(egenes.logprobs[[1]]),
          names(egenes.logprobs)))
}

#' @export
likelihoodToAdj <- function(ll.df) {
  names <- unique(as.vector(t(ll.df[,1:2])))
  r <- diag(length(names))
  dimnames(r) <- list(names, names)
  for (i in 1:nrow(ll.df)) {
    edge <- which.max(ll.df[i,3:8])
    A <- as.character(ll.df[i,1]) # as.character to force string indexing
    B <- as.character(ll.df[i,2])
    if (sum(ll.df[i,3:8] == ll.df[i,edge+2]) > 1) {
      ## multiple edges score equally well, no interaction
    } else if (edge == 1) { # A to B
      r[A,B] <- 1
    } else if (edge == 2) { # B to A
      r[B,A] <- 1
    } else if (edge == 3) { # A inhib B
      r[A,B] <- -1
    } else if (edge == 4) { # B inhib A
      r[B,A] <- -1
    } else if (edge == 5) { # A non B
      # no edges
    } else if (edge == 6) { # A eqv B
      r[A,B] <- r[B,A] <- 1
    }
  }
  return(r)
}

#' @export
scoreModel <- function(accMatrix, egenes.logprobs) {
  stopifnot(rownames(accMatrix) == names(egenes.logprobs),
            colnames(accMatrix) == names(egenes.logprobs))
  logprobs <- log.probsTo3Array(egenes.logprobs)
  # attachment likelihoods, P(D_{i,*} | M, \theta_i = j)
  # i is index over egenes (col index)
  # j is index over postive and negative attachment to S genes (row index)
  diag(accMatrix) <- 1
  upMatrix   <- accMatrix + 2 # from -1,0,1 to neg,null,pos
  downMatrix <- 2 - accMatrix # from -1,0,1 to pos,null,neg
  pickFromRow <- function(m, pick) {
    m[0:(ncol(m)-1)*nrow(m) + pick]
  }
  getParentLikelihoods <- function(distByKO, accToDistMapping) {
    # the effects are chosen by column, since we want the effects
    # as though the egene were attached to that Sgene.
    # If it was by row, then we select what other Sgenes are
    # affected in the knockout.
    apply(accToDistMapping, 2, function(KOEffect) {
      sum(pickFromRow(distByKO, KOEffect), na.rm=T)
    })
  }
  al <- apply(logprobs, 1, function(distByKO) {
    c(getParentLikelihoods(distByKO, upMatrix),
      sum(distByKO[2,], na.rm=T),
      getParentLikelihoods(distByKO, downMatrix))
  })
  rownames(al) <- c(paste("neg", rownames(accMatrix), sep="_"),
                    "unattached",
                    paste("pos", rownames(accMatrix), sep="_"))

  u <- nrow(accMatrix) + 1
  maxConnectionRatio <- apply(al[-u,], 2, max) - al[u,]
  zs <- apply(al, 2, function(x) logsum(sort(x)))
  posterior <- t(al) - zs
  ll <- sum(zs, na.rm=T)

  return(list(ll=ll, posterior=posterior, egeneLL=zs,
              maxConnectionRatio=maxConnectionRatio))
}


#' @export
scoreBestModelEstimate  <- function(egeneStruct, ..., nullRuns=0) {
  estimates <- scoreEstimates(egeneStruct, ...)
  adj <- likelihoodToAdj(estimates$scores)
  acc <- knockoutEffect(adj)
  ll <- scoreModel(acc, estimates$egenes.logprobs)
  return (c(estimates, list(adj=adj, acc=acc), ll))
}

#' @export
scrambleData <- function(e) {
  m <- matrix(sample(e$egenes), nrow=nrow(e$egenes), ncol=ncol(e$egenes),
              dimnames=dimnames(e$egenes))
  return(list(egenes=m, knockdown.cols=e$knockdown.cols, lof=e$lof))
}

#' @export
scoreBestModelEstimateNull <- function(egeneStruct, ..., nullRuns=5) {
  return(sapply(1:nullRuns, function (i) {
    scoreBestModelEstimate(scrambleData(egeneStruct), ...)$ll
  }))
}
