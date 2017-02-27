#' @export
sgenenames <- paste('s', as.vector(t(outer(c("", LETTERS),
                                           LETTERS, FUN=paste, sep=""))),
                    sep='')

#' @export
sgeneNames <- function(n) {
  if (n > length(sgenenames)) stop("too many sgenes: ", n, ", max: ",
                                   length(sgenenames))
  return(sgenenames[1:n])
}

#' @export
randomTree <- function(n, edgeGenerator=function() {1},
                       gNames=sgeneNames(n)) {
  r <- diag(n)
  rownames(r) <- colnames(r) <- gNames
  for (i in 2:n) {
    parent <- sample(i-1, size=1)
    r[parent, i] <- edgeGenerator()
  }
  r
}

#' @export
randomModelWithCycles <- function(n, edgeGenerator=function() {1},
                                  gNames=sgeneNames(n),
                                  nBackLinks=function(n) {1}) {
  if (is.function(nBackLinks)) nBackLinks <- nBackLinks(n)
  tree <- randomTree(n, edgeGenerator=edgeGenerator, gNames=gNames)
  triscalar2pair <- function(n, position) {
    if (position > n*(n-1)/2) {
      stop("triangular position is out of bounds")
    }
    i <- 1
    while(position >= n) {
      i <- i + 1
      n <- n - 1
      position <- position - n
    }
    return(c(i, position + i))
  }
  positions <- sample(n*(n-1)/2, nBackLinks)
  for (position in positions) {
    p <- triscalar2pair(n, position)
    tree[p[2],p[1]] <- edgeGenerator()
  }
  return(tree)
}
 
#' @export
posneg <- function(n=1, lower=0.75, upper=1.0, negProb=0.25) {
  sample(c(-1,1), replace=T, size=n, prob=c(negProb, 1-negProb)) *
    runif(n, min=lower, max=upper)
}

#' @export
generateEgeneParents <- function (nEgenes, sgeneNames, prior=NULL,
                                  inhibition=0.25) {
  parent <- sample(1:length(sgeneNames), size=nEgenes,
                   replace=TRUE, prob=prior)
  sign <- sample(c(-1,1), size=nEgenes, replace=TRUE,
                 prob=c(inhibition, 1 - inhibition))
  parent <- sign * parent
  parentChar <- c("n", "error", "p")[sign+2]
  names(parent) <- paste('egene', 1:nEgenes, parentChar,
                         sgeneNames[abs(parent)], sep="")
  return(parent)
}

#' @export
generateQuantization <- function(knockdown.cols, parents, acc) {
  # m vars are to account for a parent of 0, unattached
  macc <- cbind(acc, 0)
  mparents <- parents
  mparents[parents == 0] <- ncol(macc)
  sapply(knockdown.cols, function(x) {
    p <- macc[x,abs(mparents)] * -sign(parents)
    r <- (runif(parents) < abs(p)) * sign(p)
    names(r) <- names(parents)
    r
  })
}

#' @export
generateExpression <- function(quant, params) {
  r <- matrix(NA, nrow=nrow(quant), ncol=ncol(quant), dimnames=dimnames(quant))
  r[quant == -1] <- rnorm(sum(quant == -1), mean=params["neg", "mean"],
                          sd=params["neg", "sd"])
  r[quant == 0] <- rnorm(sum(quant == 0), mean=params["null", "mean"],
      sd=params["null", "sd"])
  r[quant == 1] <- rnorm(sum(quant == 1), mean=params["pos", "mean"],
      sd=params["pos", "sd"])
  return(r)
}

#' @export
generateExperiment <- function(nSgenes=8, nEgeneMult=20,
                               nEgenes=nEgeneMult*nSgenes,
                               includeSgeneEffect=FALSE,
                               nonEffects=0,
                               egeneInhibition=0.25,
                               adjGenerator=posneg,
                               modelGenerator=randomTree,
                               replicates=1,
                               params=paramGen(2, 1)) {

  if (is.function(nSgenes)) nSgenes <- nSgenes()
  if (is.function(nEgenes)) nEgenes <- nEgenes()
  adj <- modelGenerator(nSgenes, edgeGenerator=adjGenerator)
  acc <- knockoutEffect(adj)
  parents <- generateEgeneParents(nEgenes, colnames(adj),
                                  inhibition=egeneInhibition)
  if (includeSgeneEffect) {
    seffect <- 1:nrow(adj)
    names(seffect) <- rownames(adj)
    parents <- c(parents, seffect)
    nEgenes <- length(parents)
  }
  if (length(nonEffects) != 1) stop("non-scalar specified for nonEffects")
  if (nonEffects > 0) {
    ne <- rep(0, nonEffects)
    names(ne) <- paste("unattached", 1:nonEffects, sep="")
    parents <- c(parents, ne)
    nEgenes <- length(parents)
  }
  knockdown.cols <- rep(colnames(adj), each=replicates)
  quant <- generateQuantization(knockdown.cols, parents, acc)
  expr <- generateExpression(quant, params)

  r <- list(egenes=expr, knockdown.cols=knockdown.cols, lof=colnames(adj),
            sgeneAdj=adj, sgeneAcc=acc, egeneParents=parents,
            egeneQuant=quant, obsParams=params)

  r
}

#' convert an accesibility matrix into a list of the correct answers
#  1 - AtoB
#  2 - BtoA
#  3 - ArepB
#  4 - BrepA
#  5 - AnonB
#  6 - AeqvB
#  7 - A |-> B
#  8 - A <-| B
#  9 - A |-| B
#' @export
accToInteraction <- function(acc) {
  kd <- colnames(acc)
  pairs <- expand.grid(A=kd,B=kd)[as.vector(upper.tri(diag(length(kd)))),]
  pair2int <- function(p) {
    ab <- acc[p[1],p[2]]
    ba <- acc[p[2],p[1]]
    i <- (sign(ab)+1) + 3*(sign(ba)+1) + 1
    map <- c(9,  # 1 - (-1, -1) - AeqvB - 6 NA
             4,  # 2 - (0,  -1) - BrepA - 4
             8,  # 3 - (1,  -1) - AeqvB - 6 NA
             3,  # 4 - (-1,  0) - ArepB - 3
             5,  # 5 - (0,   0) - AnonB - 5
             1,  # 6 - (1,   0) - AtoB  - 1
             7,  # 7 - (-1,  1) - AeqvB - 6 NA
             2,  # 8 - (0,   1) - BtoA  - 2
             6)  # 9 - (1,   1) - AeqvB - 6
#    map <- c(NA,4,NA,3,5,1,NA,2,6)
    map[i]
#    c(max(abs(ab), abs(ba)), map[i])
  }
  r <- apply(pairs, 1, pair2int)
#  rownames(r) <- apply(pairs, 1, paste, collapse="*")
  names(r) <- apply(pairs, 1, paste, collapse="*")
  r
}

#' @export
accToLinkStrength <- function(acc) {
  kd <- colnames(acc)
  pairs <- expand.grid(A=kd,B=kd)[as.vector(upper.tri(diag(length(kd)))),]
  r <- apply(pairs, 1, function(p) {
    ab <- acc[p[1],p[2]]
    ba <- acc[p[2],p[1]]
    max(abs(ab), abs(ba))
  })
  names(r) <- apply(pairs, 1, paste, collapse="*")
  r
}

#' @export
gatherMean <- function(data) {
  n <- nrow(data$sgeneAcc)
  rep(data$obsParams[3,1], n * (n-1) / 2) 
}

#' @export
gatherFracInhibition <- function(data) {
  adj <- data$sgeneAdj
  n <- nrow(adj)
  rep(sum(adj < 0) / sum(adj != 0), n * (n-1) / 2)
}

#' @export
gatherSize <- function(data) {
  n <- nrow(data$sgeneAdj)
  rep(n, n * (n-1) / 2)
}

#' @export
egenePostToParent <- function(s) {
  p <- apply(s$posterior, 1, which.max)
  mid <- (ncol(s$posterior) + 1)/2
  sapply(p, function(i) if (i < mid) i else (mid - i))
}

#' Take in a data.frame of liklihood scores (first two columns sgene names,
#' rest are liklihoods), and make the predictions:
#'  1 - AtoB
#'  2 - BtoA
#'  3 - ArepB
#'  4 - BrepA
#'  5 - AnonB
#'  6 - AeqvB
#' @export
likelihoodToPredictions <- function(ll.df) {
  r <- apply(as.matrix(ll.df[,3:8]), 1, which.max)
  names(r) <- apply(ll.df[,1:2], 1, paste, collapse="*")
  r
}

#' @export
likelihoodToScore <- function(ll.df) {
  r <- apply(as.matrix(ll.df[,3:8]), 1, function(x) {
    -diff(sort(x, decreasing=T))[1]
    })
  names(r) <- apply(ll.df[,1:2], 1, paste, collapse="*")
  r
}

#' pairwise is the pairwise score matrix
#' egenes is the egene structure (with knockdown.cols)
#' @export
scorePairwise <- function(pairwise, egenes, params) {
  adj <- likelihoodToAdj(pairwise)
  model <- knockoutEffect(adj)
  log.probs <- exprToRegLogProbs(egenes$egenes, egenes$knockdown.cols, params)
  return(scoreModel(model, log.probs))
}

#' @export
collectRuns <- function(data, predictions, scorematrix="scores") {
  if (length(data) != length(predictions)) {
    stop("length mismatch on data and prediction lists")
  }
  tmp.sizes <- sapply(data, function(x) nrow(x$sgeneAdj))
  id <- rep(1:length(data), (tmp.sizes -1) * tmp.sizes / 2)
  size <- unlist(lapply(data, gatherSize))
  separation <- unlist(lapply(data, gatherMean))
  fracInhib <- unlist(lapply(data, gatherFracInhibition))
  strength <- unlist(lapply(data, function(d) accToLinkStrength(d$sgeneAcc)))
  model.links <- lapply(data, function(d) accToInteraction(d$sgeneAcc))
  actual <- unlist(model.links)
  pair.names <- unlist(lapply(model.links, names))
  predicted <- unlist(lapply(predictions, function(p)
                             likelihoodToPredictions(p[[scorematrix]])))
  llScore <- unlist(lapply(predictions, function(p)
                             likelihoodToScore(p[[scorematrix]])))
  correct <- actual == predicted
  TMPSTRUCT <- c(1,2,1,2,5,6,6,6,6)
  correctStructure <- TMPSTRUCT[actual] == TMPSTRUCT[predicted]

  edgeL <- c("->", "<-", "-|", "|-", ":", "<->", "|->", "<-|", "|-|")

  data.frame(id, pair.names, size, separation, fracInhib, strength,
             actual=edgeL[actual], predicted=edgeL[predicted], llScore,
             correct, correctStructure)
}

#' @export
collectEgenes <- function(data, posterior, scorematrix="scores") {
  stopifnot(length(data)==length(posterior))
  scores <- lapply(1:length(data), function(i) {
    scorePairwise(predictions[[i]][[scorematrix]],
                  data[[i]], data[[i]]$obsParams)
  })
  negenes <- sapply(scores, function(x) nrow(x$posterior))
  id <- rep(1:length(data), negenes)
  fracInhib <- rep(sapply(data, function(d) sum(d$sgeneAdj < 0) /
                          sum(d$sgeneAdj != 0)), negenes)
  ll <- rep(sapply(scores, function(x) x$ll), negenes)
  actualEgeneParent <- unlist(lapply(data, function(d) d$egeneParents))
  predEgeneParent <- unlist(lapply(scores, egenePostToParent))
  data.frame(id, fracInhib, ll, actualEgeneParent, predEgeneParent)
}

##
## Plotting
##
library(binom)

#' @export
histratio <- function(x, breaks="Sturges", conf.level=0.95, length=0.1, ...) {
  allhist <- hist(unlist(x), breaks = breaks, plot = FALSE)
  combhist <- sapply(x, function(z) hist(z, breaks = allhist$breaks,
                                           plot = FALSE)$counts)
  b <- binom.confint(combhist[,1], rowSums(combhist), methods="exact")
  x <- barplot(b[,'mean'], names=signif(allhist$mids, 2), ...)
  arrows(x, b[,'lower'], x, b[,'upper'], angle=90, code=3, length=length)
}

#' @export
multhistratio <- function(x, breaks="Sturges",
                          conf.level=0.95, length=0.1,
                          axislabelstart=NULL,
                          col=gray.colors(length(x)),
                          legend.x=NULL, legend.y=NULL, legend.txt=names(x),
                          ...) {
  allhist <- hist(unlist(x), breaks = breaks, plot = FALSE)
  combhist <- lapply(x, function(y) {
    sapply(y, function(z) hist(z, breaks = allhist$breaks,
                                           plot = FALSE)$counts)
  })
  b <- lapply(combhist, function(x) binom.confint(x[,1], rowSums(x),
                                                  methods="exact"))
  bars <- t(sapply(b, function(x) {
    a <-  x[,'mean']
    a[x[,'n']==0] <- 0
    a
  }))
  axislabel <- signif(allhist$mids, 2)
  if (!is.null(axislabelstart)) {
    axislabel[1:length(axislabelstart)] = axislabelstart
  }
  bars.x <- barplot(bars, names=axislabel,
                    beside=TRUE, col=col, ...)
  for (i in 1:length(b)) {
    lower <- b[[i]][,'lower']
    lower[b[[i]][,'n']==0] <- NA
    upper <- b[[i]][,'upper']
    upper[b[[i]][,'n']==0] <- NA
    arrows(bars.x[i,], lower, bars.x[i,], upper,
           angle=90, code=3, length=length)
  }
  if (!is.null(legend.txt) & !is.null(legend.x) & !is.null(legend.y)) {
    legend(legend.x, legend.y, names(x), fill=col)
  }
  return(b)
}
