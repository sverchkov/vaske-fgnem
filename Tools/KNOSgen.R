#!/usr/bin/env Rcmd.py
library(polynom)
library(lattice)
library(mixgauss)
library(logsum)
library(KnockoutNets)
library(KnockoutNetsSynth)
library(RCommandArgs)

usage <-"

Generate models under the KnockoutNet Synthetic model.
Options:
  OUTFILE    filename to write the R data structure of models
  MODELS     (integer) number of models to generate (default 500)
  MIN_N      (integer) minimum number of SGenes in a model (default 5)
  MAX_N      (integer) maximum number of SGenes in a model (default 15)
  EFFECTSIZE (integer) the average number of EGenes per SGene (default 20)
  S_EFFECTS  (logical) if TRUE, include Sgene effects (default FALSE)
  NONEFFECTS (integer) generate unattached genes (default 0)
  HYBREPS    (integer) number of synthetic hybridizations per knockout (def 1)
  OBS_SEP    separation between the observation peaks (default 1.5)
  OBS_SD     standard deviation of the observation peaks (default 1.0)
  BACK_FRAC  number of back links to make a cyclic model, specified as a
             fraction of N for a model (default 0.25)
  EDGE_MIN   minimum strength on an edge (default 0.75)
  INHIB      chance of generating a negative edge (default 0.25)

Example Usage:
  KNOSgen.R OUTFILE=data.RData MODELS=5 INHIB=0.75  "

outfile <- RCommandArgString("OUTFILE", errorMsg=usage)

models <- RCommandArgInteger("MODELS", default=500, gt=0, errorMsg=usage)
min.n <- RCommandArgInteger("MIN_N", default=5, gt=0, errorMsg=usage)
max.n <- RCommandArgInteger("MAX_N", default=15, gte=min.n, errorMsg=usage)
effectsize <- RCommandArgInteger("EFFECTSIZE", default=20, gte=1,
                                 errorMsg=usage)

s.effects <- RCommandArgLogical("S_EFFECTS", default=FALSE, errorMsg=usage)
noneffects <- RCommandArgInteger("NONEFFECTS", default=0, gte=0,
                                 errorMsg=usage)

hybreplicates <- RCommandArgInteger("HYBREPS", default=1, gt=0, errorMsg=usage)
mean.sep <- RCommandArgDouble("OBS_SEP", default=1.5, gt=0, errorMsg=usage)
sd <- RCommandArgDouble("OBS_SD", default=1, gt=0, errorMsg=usage)
backProportion <- RCommandArgDouble("BACK_FRAC", default=0.25, gte=0,lte=1,
                                    errorMsg=usage)
minStrength <- RCommandArgDouble("EDGE_MIN", default=0.75, gte=0, lte=1,
                                 errorMsg=usage)
inhib <- RCommandArgDouble("INHIB", default=0.25, gte=0, lte=1,
                               errorMsg=usage)
egene.inhib <- RCommandArgDouble("EFFECTINHIB", default=inhib, gte=0, lte=1,
                                 errorMsg=usage)

edgeGen <- function() posneg(negProb=inhib, lower=minStrength)

modelGen <- function(n, edgeGenerator, bp=backProportion) {
  back <- as.integer(bp*n)
  randomModelWithCycles(n, edgeGenerator=edgeGenerator, nBackLinks=back)
}

randomNSgenes <- if (min.n == max.n) {
  function() {return(max.n)}
} else {
  function() {sample(min.n:max.n, 1)}
}

data <- replicate(models,
                  generateExperiment(nSgenes=randomNSgenes,
                                     nEgeneMult=effectsize,
                                     includeSgeneEffect=s.effects,
                                     nonEffects=noneffects,
                                     egeneInhibition=egene.inhib,
                                     adjGenerator=edgeGen,
                                     modelGenerator=modelGen,
                                     params=paramGen(mean.sep, sd=sd),
                                     replicates=hybreplicates),
                  simplify=FALSE)

save(data, file=outfile)
