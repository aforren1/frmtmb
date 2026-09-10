# Lane `latent`: roxygenise frmtmb.latent in place.
source("dev/latent-env.R")
roxygen2::roxygenise("extensions/frmtmb.latent", load_code = "source")
