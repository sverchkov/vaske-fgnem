\name{RCommandArgDouble}
\alias{RCommandArgDouble}
\alias{RCommandArgInteger}
\alias{RCommandArgString}
\title{ RCommandArg Parsing }
\description{
  Reads and parses arguments from environment variables.  The argument
  can have a default argument, and numeric (double or integer) arguments
  can be optionally bounds checked.
}
\usage{
RCommandArgDouble(name, lt, lte, gt, gte, default, errorMsg)
RCommandArgInteger(name, lt, lte, gt, gte, default, errorMsg)
RCommandArgString(name, default, errorMsg)
}
\arguments{
  \item{name}{ The environment variable }
  \item{lt}{ The variable must be less than this argument }
  \item{lte}{ The variable must be less than or equal to this argument }
  \item{gt}{ The variable must be greater than this argument }
  \item{gte}{ The variable must be greater than or equal to this argument }
  \item{default}{ Default value if the environment variable is empty }
  \item{errorMsg}{ Text to be output if any checks fail }
}
\details{
  This facilitates R shell script type commands via Rcmd.py.  Arguments
  on the command line are specified with variablename=value pairs, and
  then can be read using these these functions.

  The default value is not subject to bounds checking
}
\value{
  The value in the environment variable, or the default value if the
  enviroment variable is empty.
}
\author{ Charlie Vaske <cvaske@soe.ucsc.edu> }
\examples{
  # Results in comments might correspond to the following environment
  # > RCmd.py test.R D=22.1 I=9 S="Fuzzy kittens"
  Sys.putenv(D="22.1", I="9", S="Fuzzy kittens")

  d <- RCommandArgDouble("D", lte=23, gt=0, default="11") # 22.1
  i <- RCommandArgInteger("I", lt=10, gt=8) # 9
  s <- RCommandArgString("S") # "Fuzzy Kittens"

  try(RCommandArgDouble("nonexistant")) # error
  try(RCommandArgDouble("nonexistant", default=1.0)) # 1.0

  try(RCommandArgInteger("D")) # error (D not an integer)
  try(RCommandArgDouble("S")) # error (S is not a double)
  i.asd <- RCommandArgDouble("I") # 9.0 

  try(RCommandArgDouble("D", lte=1.0)) # error (bounds)
}
\keyword{ IO }
