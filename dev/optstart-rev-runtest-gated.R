# Reviewer of lane optstart: dev/optstart-rev-runtest.R with the Stan
# tier on (the gate reads FRMTMB_BRMS_FIT_TESTS inside the process).
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-optstart/dev/stan-cache",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
source("dev/optstart-rev-runtest.R")
