\name{RCommandArgFile}
\alias{RCommandArgFile}
\title{ RCommandArg File Argument }
\description{
  Reads and opens files from an environment variable
}
\usage{
RCommandArgFile(name, default, errorMsg, open="r")
}
\arguments{
  \item{name}{ The environment variable argument name }
  \item{default}{ Default file name if the environment variable is empty }
  \item{errorMsg}{ Text to be output if any checks fail }
  \item{open}{ How to open the file }
}
\details{
  This facilitates R shell script type commands via Rcmd.py.  Arguments
  on the command line are specified with variablename=value pairs, and
  then can be read using these these functions.

  This funciton will open a file specified by the value of the
  environment variable. If there are improper permissions for reading or
  writing, or if trying to read a file that does not exist, the function
  will exit with an exception.

  See the documentation for \code{\link{file}} for \code{open} semantics.
  
}
\value{
  A connection object to the file specified by the environment
  variable.  Standard input and output can be accessed when the value of
  the environment variable is '-'.  If the value is '-' then \code{open} must
  be exactly "r" for stdin, and "w" for stdout.  If \code{open} is
  anything else, then \code{RCommandArgFile} will throw an exception.
}
\author{ Charlie Vaske <cvaske@soe.ucsc.edu> }
\examples{
  # Results in comments might correspond to the following environment
  # > RCmd.py OUT=/tmp/w IN=Makefile STDIN=- STDOUT=-
  Sys.putenv(OUT="/tmp/w", IN="../Makefile", STDIN="-", STDOUT="-")

  o <- RCommandArgFile("OUT", open="w")
  o.fn <- summary(o)$description
  close(o)
  unlink(o.fn)

  o.error <- try(RCommandArgFile("OUT", open="r"))
  stopifnot(class(o.error) == "try-error")

  i <- RCommandArgFile("IN")
  i.lines <- readLines(i)
  close(i)

  si <- RCommandArgFile("STDIN", open="r")

  so <- RCommandArgFile("STDOUT", open="w")

  std.error <- try(RCommandArgFile("STDIN", open="rw"))
  stopifnot(class(std.error) == "try-error")
  
}
\keyword{ IO }
