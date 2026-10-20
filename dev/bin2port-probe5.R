# Probe 5: once a bare multinomial() resolves its category count from
# the data (lane catmulti's change, which this lane does not depend on),
# does the adapter's multinomial path answer brms's rows log_lik:504 and
# posterior_predict:349? A stand-in for that change: brmsfamily() with a
# multinomial whose `defer` slot reads K from the response's columns, as
# the deferred categorical and simplex families read theirs.
# Log: dev/bin2port-out/probe5.txt. Seed 504.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/brmsport-bin2-adapters.R")
fake_brmsfamily <- function(family, ...) {
  if (!identical(family, "multinomial")) return(frmtmb::brmsfamily(family,
                                                                   ...))
  fam <- frmtmb::multinomial(K = 2)
  fam[["defer"]] <- function(formula, data) {
    frmtmb::multinomial(K = ncol(eval(formula[[2]], data)))
  }
  fam
}
b2 <- deparse(body(brms_adapt_cat_family))
b2 <- gsub("frmtmb::brmsfamily", "fake_brmsfamily", b2, fixed = TRUE)
body(brms_adapt_cat_family) <- parse(text = b2)[[1]]
stopifnot(any(grepl("fake_brmsfamily", deparse(body(brms_adapt_cat_family)))))
set.seed(504)
ns <- 50
nobs <- 8
ncat <- 3
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(mu1 = array(rnorm(ns * nobs), dim = c(ns, nobs)),
                   mu2 = array(rnorm(ns * nobs), dim = c(ns, nobs)))
prep$data <- list(Y = matrix(sample(1:20, nobs * ncat, TRUE), nrow = nobs),
                  ncat = ncat, trials = sample(1:20, nobs))
prep$refcat <- 1
a <- sapply(1:nobs, brms_adapt_log_lik, fun = "multinomial", prep = prep)
bp <- prep
bp$family <- brms::multinomial()
r <- sapply(1:nobs, brms:::log_lik_multinomial, prep = bp)
cat("log_lik multinomial dim", dim(a), " max abs diff vs brms",
    format(max(abs(a - r)), digits = 3), "\n")
p <- brms_adapt_posterior_predict("multinomial", 3, prep)
cat("posterior_predict multinomial dim", dim(p), " row sums == trials[3]:",
    all(rowSums(p) == prep$data$trials[3]), "\n")
