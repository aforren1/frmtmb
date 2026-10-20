# Review finding 2: the value check's comparison, brms_value_compare(),
# with the per-entry bound d <= sqrt(eps) * pmax(|ref|, 1) beside
# all.equal(). The reviewer's cases (dev/bin2port-rev-tol.R, seed 534):
# the logistic-normal block of tests.log_lik.R, one entry (the smallest
# |ref|) moved by 1e-6, 1e-4, 1e-2 and 1, and the catlpdf mutation of
# dev/bin2port-rev-mutate.R (every log density + 1e-6), which the old
# rule let through at log_lik:534. A control: the unmoved values hold.
# Log: dev/bin2port-out/tol.txt.
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
old_rule <- function(cur, ref) isTRUE(all.equal(ref, cur))
show <- function(label, x) {
  r <- brms_value_compare(x, ref)
  cat(sprintf("%-44s all.equal alone: %-5s  new rule: %-5s  %s\n", label,
              old_rule(x, ref), r$held, r$msg))
}
cat("|ref| range:", format(range(abs(ref)), digits = 3), "\n")
show("control: the adapter's values", cur)
k <- which.min(abs(ref))
for (eps in c(1e-6, 1e-4, 1e-2, 1)) {
  bad <- cur
  bad[k] <- bad[k] + eps
  show(sprintf("one entry (|ref| %.3g) + %g", abs(ref[k]), eps), bad)
}
show("every entry + 1e-6 (the catlpdf mutation)", cur + 1e-6)
