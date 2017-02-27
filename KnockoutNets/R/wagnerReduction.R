##
## WARNING!  this is as naive, fragile algorithm which
## recurses infinitely on cycles. (big big hack)
##
wagnerReduction <- function(accMatrix) {
  stopifnot(is.numeric(accMatrix))
  stopifnot(length(dim(accMatrix)) == 2)
  stopifnot(nrow(accMatrix) == ncol(accMatrix))

  acc <- accMatrix != 0
  diag(acc) <- FALSE
  adj <- acc

  visited <- logical(nrow(acc))

  pruneAcc <- function(i) {
    for (j in which(acc[i,])) {
      if (!any(acc[j,])) {
        visited[j] <<- TRUE
      } else {
        pruneAcc(j)
      }
    }
    for (j in which(acc[i,])) {
      adj[i,adj[j,] & acc[i,]] <<- FALSE
    }
    visited[i] <<- TRUE
  }
  for (i in 1:nrow(adj)) {
    if (!visited[i]) pruneAcc(i)
  }
  diag(adj) <- TRUE
  accMatrix[!adj] <- 0
  return(accMatrix)
}


identicalColumns <- function(df) {
  return(match(df, df))
}

mergeEquiv <- function(accMatrix, sep=",") {
  stopifnot(is.numeric(accMatrix))
  stopifnot(length(dim(accMatrix)) == 2)
  stopifnot(nrow(accMatrix) == ncol(accMatrix))
  if (is.null(dimnames(accMatrix))) {
    l <- paste("g", 1:nrow(accMatrix), sep='')
    dimnames(accMatrix) <- list(l,l)
  }

  map <- identicalColumns(data.frame(accMatrix))
  n <- sapply(split(1:length(map), map), function(x) {
    paste(colnames(accMatrix)[x], collapse="|")
  })
  umap <- unique(map)
  acc <- accMatrix[umap,umap, drop=F]
  dimnames(acc) <- list(n[as.character(umap)], n[as.character(umap)])
  return(acc)
}

accToDot <- function(accMatrix, forwardLinksOnly=TRUE) {
  stopifnot(is.numeric(accMatrix))
  stopifnot(length(dim(accMatrix)) == 2)
  stopifnot(nrow(accMatrix) == ncol(accMatrix))
  if (is.null(dimnames(accMatrix))) {
    l <- paste("g", 1:nrow(accMatrix), sep='')
    dimnames(accMatrix) <- list(l,l)
  }

  nodes <- paste("   \"",rownames(accMatrix),"\";", sep="")
  links <- NULL
  if (nrow(accMatrix) > 1) {
    acc <- sign(accMatrix)
    
    kd <- colnames(acc)
    pairs <- expand.grid(A=kd,B=kd)[as.vector(upper.tri(diag(length(kd)))),]
    
    r <- t(apply(pairs, 1, function(p) c(acc[p[1],p[2]],acc[p[2],p[1]])))
    
    nonint <- apply(abs(r), 1, sum) == 0
    if (any(!nonint)) {
      pairs <- pairs[!nonint,,drop=FALSE]
      r <- r[!nonint,,drop=FALSE]
      
      if (forwardLinksOnly) {
        swapAB <- apply(r, 1, function(x) (x[1] == 0 & x[2] != 0))
        pairs[swapAB,] <- pairs[swapAB,2:1]
        r[swapAB,] <- r[swapAB,2:1]
      }
      
      arrowspec <- c("tee", "none", "normal")
      
      links <- paste("   \"", pairs[,1], "\" -> \"", pairs[,2],
                     "\" [arrowtail=", arrowspec[r[,2]+2],
                     ",arrowhead=", arrowspec[r[,1]+2],
                     "] ;", sep="")
    }
  }
  return(c("digraph {", nodes, links, "}"))
}
