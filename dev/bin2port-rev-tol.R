# Reviewer: how large an error the value check's all.equal() lets
# through on the logistic-normal block of tests.log_lik.R (log_lik:534),
# whose values span many orders of magnitude. Seed 534.
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
cat("|ref| quantiles:", format(quantile(abs(ref), c(0, .1, .5, .9, 1)),
                               digits = 3), "\n")
for (eps in c(1e-6, 1e-4, 1e-2, 1)) {
  bad <- cur
  k <- which.min(abs(ref))
  bad[k] <- bad[k] + eps
  cat(sprintf(paste("one entry (|ref| = %.3g) off by %g: all.equal holds",
                    "%s; max elementwise rel %.3g\n"), abs(ref[k]), eps,
              isTRUE(all.equal(ref, bad)),
              max(abs(bad - ref) / pmax(abs(ref), 1e-300))))
}
