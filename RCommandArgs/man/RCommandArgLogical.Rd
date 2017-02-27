\name{RCommandArgLogical}
\alias{RCommandArgLogical}
\title{ RCommandArg Logical Argument }
\description{
  Reads logical (i.e. boolean/binary) variables from an environment variable
}
\usage{
RCommandArgLogical(name, default, errorMsg)
}
\arguments{
  \item{name}{ The environment variable argument name }
  \item{default}{ Default value if the environment variable is empty }
  \item{errorMsg}{ Text to be output if value is absent }
}
\details{
  Reads a logical (Boolean) variable from an environment variable.
}
\value{
  Returns \code{TRUE} or \code{FALSE}, or an error if the string can not
  be coerced to a logical value.  If the string is not an R logical
  value, coercion to a double value and then to logical is attempted.
}
\author{ Charlie Vaske <cvaske@soe.ucsc.edu> }
\examples{
  # Values that should work
  # > RCmd.py ST=T SF=F LTRUE=TRUE LFALSE=FALSE E0=0 E1=1
  Sys.putenv(ST="T", SF="F", LTRUE="TRUE", LFALSE="FALSE", E0="0", E1="1")
  stopifnot(RCommandArgLogical("ST") == TRUE)
  stopifnot(RCommandArgLogical("LTRUE") == TRUE)
  stopifnot(RCommandArgLogical("E1") == TRUE)
  stopifnot(RCommandArgLogical("SF") == FALSE)
  stopifnot(RCommandArgLogical("LFALSE") == FALSE)
  stopifnot(RCommandArgLogical("E0") == FALSE)

  # These should not work
  Sys.putenv(GIB="aoeu", UCGIB="AOEU")
  stopifnot(class(try(RCommandArgLogical("GIB"))) == "try-error")
  stopifnot(class(try(RCommandArgLogical("UCGIB"))) == "try-error")
}
\keyword{ IO }
