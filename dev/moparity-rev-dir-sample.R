# Reviewer: frm_sample() with a simo prior below 1/2 on the formula
# route (frm() refuses it for a MAP fit), and on the fit route of a MAP
# fit; as_tmbstan()'s message on a MAP fit with and without a simo row.
# Usage: Rscript dev/moparity-rev-dir-sample.R
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(3)
n <- 120
dm <- data.frame(inc = sample(0:3, n, TRUE), z = rnorm(n))
dm$y <- 1 + 0.6 * c(0, 1, 1.4, 1.5)[dm$inc + 1] + 0.3 * dm$z + rnorm(n)
pr04 <- set_prior("dirichlet(0.4)", class = "simo", coef = "moinc1")
r <- tryCatch(frm_sample(bf(y ~ mo(inc) + z), family = gaussian(), data = dm,
                         prior = pr04, chains = 2, iter = 600, seed = 1,
                         refresh = 0),
              error = function(e) e)
cat("formula route, dirichlet(0.4):",
    if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else {
      M <- r$draws
      paste("ok; mean simo", paste(format(colMeans(M[, grep("^simo", colnames(M))]),
                                          digits = 3), collapse = " "))
    }, "\n")
pr2 <- set_prior("dirichlet(2)", class = "simo", coef = "moinc1")
fit <- frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm, prior = pr2)
r2 <- tryCatch(suppressWarnings(frm_sample(fit, chains = 2, iter = 600,
                                           seed = 1, refresh = 0)),
               error = function(e) e)
cat("fit route of a dirichlet(2) MAP fit:",
    if (inherits(r2, "error")) paste("ERROR", conditionMessage(r2)) else {
      M <- r2$draws
      paste("ok; mean simo", paste(format(colMeans(M[, grep("^simo", colnames(M))]),
                                          digits = 3), collapse = " "))
    }, "\n")
msgs <- character(0)
r3 <- withCallingHandlers(tryCatch(as_tmbstan(fit, chains = 1, iter = 200,
                                              refresh = 0, seed = 1),
                                   error = function(e) e),
                          message = function(m) {
                            msgs <<- c(msgs, conditionMessage(m))
                            invokeRestart("muffleMessage")
                          }, warning = function(w) invokeRestart("muffleWarning"))
cat("as_tmbstan on MAP with simo row: messages", length(msgs), ":",
    paste(substr(msgs, 1, 100), collapse = " | "), "\n")
fit0 <- frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm,
            prior = set_prior("normal(0, 5)", class = "b"))
msgs <- character(0)
r4 <- withCallingHandlers(tryCatch(as_tmbstan(fit0, chains = 1, iter = 200,
                                              refresh = 0, seed = 1),
                                   error = function(e) e),
                          message = function(m) {
                            msgs <<- c(msgs, conditionMessage(m))
                            invokeRestart("muffleMessage")
                          }, warning = function(w) invokeRestart("muffleWarning"))
cat("as_tmbstan on MAP without simo row: messages", length(msgs), ":",
    paste(substr(msgs, 1, 100), collapse = " | "), "\n")
