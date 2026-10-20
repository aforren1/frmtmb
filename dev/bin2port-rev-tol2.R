# Reviewer, fix round 1: the cases of dev/bin2port-rev-tol.R put through
# the lane's NEW brms_value_compare() (dev/brmsport-bin2-adapters.R).
# Seed 534. Log dev/bin2port-rev-out/tol2.txt.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/brmsport-bin2-adapters.R")
set.seed(534)
ns <- 50; nobs <- 8; ncat <- 3
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
Y <- matrix(sample(1:20, nobs * ncat, TRUE), nrow = nobs)
prep$data <- list(Y = Y / rowSums(Y), ncat = ncat)
prep$refcat <- 1
prep$dpars <- list(mu2 = rnorm(ns), mu3 = rnorm(ns),
                   sigma2 = rexp(ns, 10), sigma3 = rexp(ns, 10))
bprep <- prep
bprep$family <- brms::brmsfamily("logistic_normal")
ref <- sapply(1:nobs, brms:::log_lik_logistic_normal, prep = bprep)
cur <- sapply(1:nobs, brms_adapt_log_lik, fun = "logistic_normal",
              prep = prep)
show <- function(label, x) {
  r <- brms_value_compare(x, ref)
  cat(sprintf("%-44s held %-5s  %s\n", label, r$held, r$msg))
}
show("unmoved", cur)
k <- which.min(abs(ref))
for (eps in c(1e-6, 1e-4, 1e-2, 1)) {
  bad <- cur; bad[k] <- bad[k] + eps
  show(sprintf("entry |ref| = %.3g moved by %g", abs(ref[k]), eps), bad)
}
show("every entry + 1e-6 (catlpdf)", cur + 1e-6)
kb <- which.max(abs(ref))
bad <- cur; bad[kb] <- bad[kb] * (1 + 1e-7)
show(sprintf("largest entry (%.3g) x (1 + 1e-7)", ref[kb]), bad)
bad <- cur; bad[kb] <- bad[kb] * (1 + 1e-9)
show(sprintf("largest entry (%.3g) x (1 + 1e-9)", ref[kb]), bad)
