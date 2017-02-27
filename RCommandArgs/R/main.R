usage <- function(..., usage) {
  cat(usage, "\n\n")
  stop(...)
}

"RCommandArgDouble" <-
function(name, lt, lte, gt, gte, default, errorMsg) {
  if (missing(errorMsg)) errorMsg <- ''

  stringval <- Sys.getenv(name)
  if (stringval == "") {
    if (missing(default)) {
      usage("Parameter ", name, " unset, with no default. ", usage=errorMsg)
    } else {
      return(default)
    }
  }
  doubleval <- as.double(stringval)
  if (is.na(doubleval)) {
    usage("Argument ", name, " must be a double. ", usage=errorMsg)
  }

  if (!missing(lt) && !(doubleval < lt)) {
    usage("Argument ", name, " must be less than ", lt, ". ", usage=errorMsg)
  }
  if (!missing(lte) && !(doubleval <= lte)) {
    usage("Argument ", name, " must be less than or equal to ", lte, ". ",
          usage=errorMsg)
  }
  if (!missing(gt) && !(doubleval > gt)) {
    usage("Argument ", name, " must be greater than ", gt, ". ",
          usage=errorMsg)
  }
  if (!missing(gte) && !(doubleval >= gte)) {
    usage("Argument ", name, " must be greater than or equal to ",
         gte, ". ", usage=errorMsg)
  }
  return(doubleval)
}

"RCommandArgInteger" <-
function(name, lt, lte, gt, gte, default, errorMsg) {
  if (missing(errorMsg)) errorMsg <- ''

  stringval <- Sys.getenv(name)
  if (stringval == "") {
    if (missing(default)) {
      usage("Parameter ", name, " unset, with no default. ",
            usage=errorMsg)
    } else {
      return(default)
    }
  }
  integerval <- as.integer(stringval)
  if (stringval != as.character(integerval)) {
    usage("Argument ", name, " must be a integer. ",
          usage=errorMsg)
  }

  if (!missing(lt) && !(integerval < lt)) {
    usage("Argument ", name, " must be less than ", lt, ". ",
          usage=errorMsg)
  }
  if (!missing(lte) && !(integerval <= lte)) {
    usage("Argument ", name, " must be less than or equal to ", lte, ". ",
          usage=errorMsg)
  }
  if (!missing(gt) && !(integerval > gt)) {
    usage("Argument ", name, " must be greater than ", gt, ". ",
          usage=errorMsg)
  }
  if (!missing(gte) && !(integerval >= gte)) {
    usage("Argument ", name, " must be greater than or equal to ",
         gte, ". ", usage=errorMsg)
  }
  return(integerval)
}

"RCommandArgLogical" <-
function(name, default, errorMsg) {
  if (missing(errorMsg)) errorMsg <- ''

  stringval <- Sys.getenv(name)
  if (stringval == "") {
    if (missing(default)) {
      usage("Missing required parameter ", name, ". ",
            usage=errorMsg)
    } else {
      return(default)
    }
  }
  lval <- as.logical(stringval)
  if (is.na(lval)) { # attempt to switch from number to logical
    dval <- as.double(stringval)
    if (is.na(dval)) {
      usage("Parameter ", name, " with value ", stringval,
           " can't be coerced to logical. ",
            usage=errorMsg)
    } else {
      return(as.logical(dval))
    }
  } else {
    return(lval)
  }
}

"RCommandArgString" <-
function(name, default, errorMsg) {
  if (missing(errorMsg)) errorMsg <- ''

  stringval <- Sys.getenv(name)
  if (stringval == "") {
    if (missing(default)) {
      usage("Missing parameter ", name, " with no default. ", usage=errorMsg)
    } else {
      return(default)
    }
  }
  return(stringval)
}

"RCommandArgSwitch" <-
function(name, default, table, errorMsg) {
  if (missing(errorMsg)) errorMsg <- ''

  stringval <- Sys.getenv(name)
  if (stringval == "") {
    if (missing(default)) {
      usage("Missing parameter ", name, " with no default. ", usage=errorMsg)
    } else {
      stringval <- default
    }
  }
  stringval <- match.arg(stringval, names(table))
  return(table[[stringval]])
}

"RCommandArgFile" <-
function(name, default, errorMsg, open="r") {
  if (missing(errorMsg)) errorMsg <- ''

  stringval <- Sys.getenv(name)
  if (stringval == "") {
    if (missing(default)) {
      usage("Missing parameter ", name, " with no deafult. ", usage=errorMsg)
    } else {
      stringval <- default
    }
  }

  if (stringval == "-") {
    if (open =="r") {
      return(stdin())
    } else if (open == "w") {
      return(stdout())
    } else {
      usage("open must be either 'r' or 'w' when file is '-' to specify ",
           "either stdin or stdout.", usage=errorMsg)
    }
  }

  fileval <- file(description=stringval, open=open)
  return(fileval)
}
