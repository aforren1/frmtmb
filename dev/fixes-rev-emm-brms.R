# Reviewer of lane fixes, claim 1: brms 2.23.0's emmeans() at frmtmb's
# estimates on cases the lane did not run: ns(), a log() in an
# interaction with a factor under at =, and a transform in sigma.
#   Rscript dev/fixes-rev-emm-brms.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/fixes-rev-stan-cache"))
suppressMessages({library(testthat); library(frmtmb); library(emmeans)
  library(splines)})
cat("LIB", find.package("frmtmb"), "\n")
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = globalenv())
  brms_rename <- frmtmb:::brms_rename
}
set.seed(7101)
n <- 150
d <- data.frame(x = rnorm(n), z = runif(n, 0.5, 4), w = runif(n, 1, 3),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$yc <- rpois(n, d$w * exp(0.2 + 0.3 * log(d$z) + 0.2 * (d$f == "b")))
d$yg <- 1 + 0.5 * d$z - 0.1 * d$z^2 + as.numeric(d$f) +
  0.3 * log(d$z) * (d$f == "c") + rnorm(n, 0, exp(-0.3 + 0.2 * d$z))
fmt <- function(v) paste(format(v, digits = 9), collapse = " ")
run <- function(lab, bform, fam, specs, ...) {
  cat("==", lab, "\n")
  fit <- frm(bform, data = d, family = fam)
  bb <- tryCatch(brms_fixed_fit(bform, fam, d, fit, ndraws = 4),
                 error = function(e) conditionMessage(e))
  if (is.character(bb)) { cat("  brms fixed fit ERROR", bb, "\n"); return() }
  a <- tryCatch(summary(emmeans(fit, specs, ...))$emmean,
                error = function(e) conditionMessage(e))
  b <- tryCatch(summary(emmeans(bb, specs, ...))$emmean,
                error = function(e) conditionMessage(e))
  cat(sprintf("  %-7s %s\n  %-7s %s\n", "frmtmb", fmt(a), "brms", fmt(b)))
  if (is.numeric(a) && is.numeric(b)) {
    cat(sprintf("  max rel diff %.3g\n", max(abs(a - b) / abs(b))))
  }
}
run("ns(z, 3) + f at z = c(1, 5)", bf(yc ~ ns(z, 3) + f), poisson(),
    ~ z + f, at = list(z = c(1, 5)))
run("log(z) * f | z at z = c(1, 2)", bf(yg ~ log(z) * f), gaussian(),
    ~ f | z, at = list(z = c(1, 2)))
run("sigma ~ scale(z), dpar = sigma, at z = c(1, 3)",
    bf(yg ~ f, sigma ~ scale(z)), gaussian(), ~ z, dpar = "sigma",
    at = list(z = c(1, 3)))
run("poly(z, 2) + f + offset(log(w)) at z = 6", bf(yc ~ poly(z, 2) + f +
                                                    offset(log(w))),
    poisson(), ~ f, at = list(z = 6))
