#!/usr/bin/python
"""\
Rcmd.py - wrapper for calling an R script from the command line, with
named arguments being turned into environment variables.  Use the
RCommandArgs library in R to retrieve these arguments.

Example Script:

   #!/usr/bin/env Rcmd.py
   library(RCommandArgs)
   infile <- RCommandArgString(\"IN\")
   cutoff <- RCommandArgDouble(\"CUTOFF\")
   ...

Example usage:

  example.R IN=data.tab CUTOFF=0.5
"""

import os, sys, re, stat

usage = __doc__

if (len(sys.argv) < 2):
    print usage
    sys.exit(1)
    
rscriptfile = sys.argv[1]
rscript = ""

if os.access(rscriptfile, os.R_OK):
    rscript = rscriptfile
else:
    print "Command script " + rscriptfile + " is not readable"
    sys.exit(2)

def shellquote(s):
    return "\"" + re.sub("\"", "\\\"", re.sub("\$", "\$", s)) + "\""

args = " ".join([shellquote(a) for a in sys.argv[2:] if re.match("^-",a)])
env = " ".join([shellquote(a) for a in sys.argv[2:] if not re.match("^-",a)])
cmd =  "env " + env + " R --vanilla --slave " + args + " --file=" + rscript
os.system(cmd)
