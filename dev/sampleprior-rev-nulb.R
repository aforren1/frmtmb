# Reviewer of lane sampleprior: an identity-link nu with brms's lb = 1
# given on the prior: the bound reaches tmbstan.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
fam <- frmtmb::brmsfamily("student", link_nu = "identity")
od <- cap_obj(y ~ x, D$real, fam,
              prior = set_prior("gamma(2, 0.1)", class = "nu", lb = 1))
cat("par:", names(od$par), "\nlower:", CAP$args$lower, "\nupper:", CAP$args$upper, "\n")
od <- cap_obj(y ~ x, D$real, fam)
cat("default: lower", if (is.null(CAP$args$lower)) "NULL (no bound)" else CAP$args$lower, "\n")
