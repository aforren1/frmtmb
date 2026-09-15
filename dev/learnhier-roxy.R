# Lane `learnhier`: roxygenise frmtmb.learn in place.
source("dev/learnhier-env.R")
roxygen2::roxygenise("extensions/frmtmb.learn", load_code = "source")
