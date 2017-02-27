\name{RCommandArgs-package}
\alias{RCommandArgs-package}
\alias{RCommandArgs}
\docType{package}
\title{
RCommandArgs Enviroment Variable Argumens
}
\description{
Reads arguments from environment variables
}
\details{
\tabular{ll}{
Package: \tab RCommandArgs\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-07-15\cr
License: \tab Suspended\cr
}

This package reads an environment variable, parses it into the
appropriate type, and checks bounds.  The user can supply error
messages, and default values.

}
\author{
Charlie Vaske <cvaske@soe.ucsc.edu>
}
\keyword{ package }

\seealso{

}
\examples{

  Sys.putenv(PRIORPROB=0.5)
  d <- RCommandArgDouble("PRIORPROB", lte=1.0, gte=0.0)

  Sys.putenv(MARGINOP="sum")
  marginop <- RCommandArgSwitch("MARGINOP", default=max,
                                table=list(max=max, sum=sum, avg=mean))
}
