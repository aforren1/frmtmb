# Lane fixes, punch round, m8: what conditional_effects(method =
# "posterior_predict") reports as estimate__ on models without
# truncation, brms 2.23.0 at frmtmb's estimates (Fixed_param draws)
# beside frmtmb: brms's is the median of its predictive draws
# (robust = TRUE), frmtmb's the expected response.
#   Rscript dev/fixes-ce-pred.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
suppressMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
henv <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = henv)
}
set.seed(808)
n <- 200
d <- data.frame(x = runif(n, -1, 1))
d$yg <- 1 + d$x + rnorm(n)
d$yp <- rpois(n, exp(1 + 0.5 * d$x))
d$yl <- exp(0.5 + 0.4 * d$x + rnorm(n, 0, 0.6))
xs <- c(-1, 0, 1)
fmt <- function(v) paste(format(v, digits = 6), collapse = " ")
for (cs in list(list(yg ~ x, gaussian(), stats::gaussian()),
                list(yp ~ x, poisson(), stats::poisson()),
                list(yl ~ x, lognormal(), brms::lognormal()))) {
  fit <- frm(bf(cs[[1]]), data = d, family = cs[[2]])
  bb <- henv$brms_fixed_fit(brms::bf(cs[[1]]), cs[[3]], d, fit,
                            ndraws = 4000)
  set.seed(1)
  b <- brms::conditional_effects(bb, "x", method = "posterior_predict",
                                 int_conditions = list(x = xs))[[1]]
  f <- suppressMessages(conditional_effects(
    fit, "x", method = "posterior_predict", ndraws = 4000,
    int_conditions = list(x = xs)))[[1]]
  e <- suppressMessages(conditional_effects(
    fit, "x", method = "posterior_epred", int_conditions = list(x = xs)))[[1]]
  cat("==", deparse(cs[[1]]), cs[[2]]$family, "\n")
  cat("  brms predict estimate__ :", fmt(b$estimate__), "\n")
  cat("  frm  predict estimate__ :", fmt(f$estimate__), "\n")
  cat("  frm  epred estimate__   :", fmt(e$estimate__), "\n")
  cat("  brms predict lower/upper:", fmt(b$lower__), "|", fmt(b$upper__), "\n")
  cat("  frm  predict lower/upper:", fmt(f$lower__), "|", fmt(f$upper__), "\n")
}
