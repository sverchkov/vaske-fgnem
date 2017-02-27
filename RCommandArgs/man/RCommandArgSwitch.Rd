\name{RCommandArgSwitch}
\title{ RCommandArg Switch Argument }
\description{
  Reads a string as an argument, and performs lookup in a table to get
  an R value
}
\usage{
RCommandArgSwitch(name, default, table, errorMsg)
}
\arguments{
  \item{name}{ The environment variable argument name }
  \item{default}{ Default file name if the environment variable is empty }
  \item{table}{ Lookup table for the argument }
  \item{errorMsg}{ Text to be output if any checks fail }
}
\details{
  This facilitates R shell script type commands via Rcmd.py.  Arguments
  on the command line are specified with variablename=value pairs, and
  then can be read using these these functions.

  This function will return pre-defined R values for strings.  This is
  useful, for example, if you want to specify different pre-defined
  functions in a script. The table is simply a list object, and the
  list element's names are the keys of the table.
  
}
\value{
  The R value specified in the table.  If the environment variable or
  default value are not keys in the table, execution stops.
}
\author{ Charlie Vaske <cvaske@soe.ucsc.edu> }
\examples{
  # Results in comments might correspond to the following environment
  # > RCmd.py test.R D=22.1 I=9 S="Fuzzy kittens"
  Sys.putenv(MARGINOP="avg")

  marginop <- RCommandArgSwitch("MARGINOP", default="max",
                                table=list(max=max, sum=sum, avg=mean))

  marginop(1:4) # should be 2.5

  Sys.putenv(MARGINOP="aoeu")  
  try(marginop <- RCommandArgSwitch("MARGINOP", default="max",
                                table=list(max=max, sum=sum, avg=mean)))
  # should fail, since aoeu isn't in the table
}
\keyword{ IO }
