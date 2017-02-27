#' @export
isValidScoreMatrix <- function(matrix) {
  if (length(dim(matrix)) !=2 || ncol(matrix) != 8) {
    return(false)
  }
  n <- length(idsOfScoreMatrix(matrix))
  return(nrow(matrix) == (n*(n-1)/2))
}

#' @export
idsOfScoreMatrix <- function(scoreDataframe) {
  # the transpose is necessary to ensure that a < b always
  return(unique(as.vector(t(as.matrix(scoreDataframe)[,1:2]))))
}


##' Returns a matrix of indices, such that the first two indices
##' specify the edge pair, the triple identifies the transitivity
##' factor.
##'
##' rather ungraceful :(
##' @export
edgeAssignments <- function(n) {
  result <- matrix(NA, nrow=(n*(n-1)*(n-2)/2), ncol=3)
  idx <- 1
  for (i in 1:(n-2)) {
    for (j in (i+1):(n-1)) {
      for (k in (j+1):n) {
        if (i != k && j != k) {
          result[idx + 0,] <- c(i, j, k)
          result[idx + 1,] <- c(i, k, j)
          result[idx + 2,] <- c(j, k, i)
          idx <- idx + 3
        }
      }
    }
  }
  return(result)
}

##' Take in a score matrix, and return a matrix, where each column
##' consists of the neigboring message's row indices in the message
##' matrix
##' @export
variableNeighbors <- function(scoreDF, edgeAssignments) {
  stringIds <- idsOfScoreMatrix(scoreDF)
  ids <- matrix(match(as.matrix(scoreDF)[,1:2], stringIds), ncol=2)
  result <-   apply(ids, 1, function(x) {
    which(edgeAssignments[,1] == x[1] & edgeAssignments[,2] == x[2])
  })
  if (is.null(dim(result))) {
    dim(result) <- c(1,length(result))
  }
  return(result)
}

#' @export
messagesToFactor <- function(constants, msgsToVars, varNeighbors, norm=NULL) {
  if (nrow(msgsToVars) == 0 || ncol(msgsToVars) != ncol(constants)) {
    stop("malformed msgsToVars")
  }
  if (ncol(varNeighbors) != nrow(constants)) {
    stop("malformed varNeighbors")
  }

  result <- matrix(NA, nrow=nrow(msgsToVars), ncol=ncol(msgsToVars))

  for (i in 1:nrow(constants)) {
    outgoing <- varNeighbors[,i, drop=FALSE]
    total <- constants[i,] + colSums(msgsToVars[outgoing,,drop=FALSE])
    for (msg in outgoing) {
      result[msg,] <- total - msgsToVars[msg,,drop=FALSE]
    }
  }
  if (!is.null(norm)) {
    result <- result - apply(result,1,norm)
  }
  return(result)
}

#' @export
totalVarBelief <- function(constants, msgsToVars, varNeighbors) {
  result <- matrix(NA, nrow=nrow(constants), ncol=ncol(constants))
  for (i in 1:nrow(constants)) {
    result[i,] <- constants[i,,drop=FALSE] +
      colSums(msgsToVars[varNeighbors[,i,drop=FALSE],,drop=FALSE])
  }
  return(result)
}

##' Return a 3D array, with axes B:C, A:C, and A:B, where each axis
##' takes on one of the 6 edges types
##' @export
generateTransitivityPotential <- function() {
  ## order is -> <- -| |- : <->
  n <- c("->", "<-", "-|", "|-", ":", "<->")
  f <- c(1,0,-1, 0,0,1)
  b <- c(0,1, 0,-1,0,1)
  
  result <- matrix(rep(NA, 6*6*6))
  dim(result) <- c(6,6,6)
  dimnames(result) <- list(paste("B", n, "C", sep=''),
                           paste("A", n, "C", sep=''),
                           paste("A", n, "B", sep=''))
  
  for(i in 1:6) {
    Mab <- f[i]
    Mba <- b[i]
    for(j in 1:6) {
      Mac <- f[j]
      Mca <- b[j]
      for (k in 1:6) {
        Mbc <- f[k]
        Mcb <- b[k]
        notTrans <-function(ab, bc, ac) {
          return (ab & bc & (ac != ab*bc))
        }
        chains <- matrix(c(Mab, Mbc, Mac, # 6 possible permutations
                           Mac, Mcb, Mab,
                           Mba, Mac, Mbc,
                           Mbc, Mca, Mba,
                           Mca, Mab, Mcb,
                           Mcb, Mba, Mca), nrow=3)
        mismatch <- any(apply(chains, 2, function(x) notTrans(x[1],x[2],x[3])))
        result[i,j,k] <- if (mismatch) -Inf else 0
      }
    }
  }
  return(result)
}

#' @export
cross <- function(a, b) {
  return(as.vector(outer(b, a, FUN="+")))
}

#source("~/src/cvaske/R/logsum.R")

#' @export
messagesToVars <- function(potential, msgsToFactors, norm=T) {
  if (nrow(msgsToFactors) %% 3 != 0 || ncol(msgsToFactors) != 6) {
    stop("malformed msgsToFactors")
  }
  result <- matrix(NA, nrow=nrow(msgsToFactors), ncol=ncol(msgsToFactors))
  for (i in seq(1, nrow(msgsToFactors), by=3)) {
    ab <- msgsToFactors[i,]
    ac <- msgsToFactors[i+1,]
    bc <- msgsToFactors[i+2,]

    bcOac <- outer(bc, ac, FUN="+")
    result[i,] <- apply(potential, 3, function(x) {max(x + bcOac)})

    bcOab <- outer(bc, ab, FUN="+")
    result[i+1,] <- apply(potential, 2, function(x) {max(x + bcOab)})

    acOab <- outer(ac, ab, FUN="+")
    result[i+2,] <- apply(potential, 1, function(x) {max(x + acOab)})
  }
  if (norm) {
    weights <- apply(result,1,logsum)
  } else {
    weights <- apply(result,1,median)
  }
  result <- result - weights
  return(result)
}

#' @export
passMessages <- function(scoreDF, damp=0.5, iters=5, normToV=T, normToF=NULL) {
  if (!isValidScoreMatrix(scoreDF)) {
    stop("invalid score matrix")
  }
  if (nrow(scoreDF) < 3) {
    return(list(scores=scoreDF))
  }
  ids <- idsOfScoreMatrix(scoreDF)
  ea <- edgeAssignments(length(ids))
  vn <- variableNeighbors(scoreDF, ea)
  constants <- as.matrix(scoreDF[,3:ncol(scoreDF)])
  transFac <- generateTransitivityPotential()

  msgsToVars <- matrix(0, nrow=nrow(ea), ncol=6)
  msgsToFactors <- matrix(0, nrow=nrow(ea), ncol=6)

  oldToV <- msgsToVars
  oldToF <- msgsToFactors
  
  fDelta <- rep(0, iters)
  vDelta <- rep(0, iters)
  fDelta2 <- rep(0, iters)
  vDelta2 <- rep(0, iters)
  for (i in 1:iters) {
    newToF <- (1-damp)*msgsToFactors + damp * messagesToFactor(constants,
                                                               msgsToVars, vn,
                                                               norm=normToF)
    fDelta[i] <- sum((newToF - msgsToFactors)^2)
    fDelta2[i] <- sum((newToF - oldToF)^2)
    oldToF <- msgsToFactors
    msgsToFactors <- newToF

    newToV <- (1-damp)*msgsToVars + damp *
      messagesToVars(transFac, msgsToFactors, norm=normToV)
    vDelta[i] <- sum((newToV - msgsToVars)^2)
    vDelta2[i] <- sum((newToV - oldToF)^2)
    oldToV <- msgsToVars
    msgsToVars <- newToV
  }

  final <- cbind(scoreDF[,1:2], totalVarBelief(constants, msgsToVars, vn))
  colnames(final) <- colnames(scoreDF)

  result <- list(scores=final,
                 toVars=msgsToVars, toFactors=msgsToFactors,
                 fDelta=fDelta, vDelta=vDelta,
                 fDelta2=fDelta2, vDelta2=vDelta2)
}
