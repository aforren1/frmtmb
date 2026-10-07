# How many warnings does posterior_samples() on draws raise, and does
# expect_warning() let one escape? (frmtmb.sample
# test-brms-shapes-draws.R:106 escaped one on the Ubuntu runner.)
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb.sample); library(testthat)})
set.seed(1)
dd <- data.frame(x = rnorm(30)); dd$y <- rnorm(30, 1 + dd$x)
fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
M <- matrix(rep(est, each = 6) + rnorm(6 * length(est), 0, 0.05), 6,
            dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
ds <- structure(list(stanfit = NULL,
                     draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                                   lp__ = 0), fit = fit),
                class = "frmtmb_draws")
for (who in c("posterior_samples", "brms::posterior_samples")) {
  f <- eval(parse(text = who))
  w <- character()
  withCallingHandlers(f(ds), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  cat(who, ": ", length(w), " warnings; generic from ",
      environmentName(environment(f)), "\n", sep = "")
}
local_edition(3)
esc <- 0
withCallingHandlers(
  expect_warning(brms::posterior_samples(ds), "deprecated"),
  warning = function(x) { esc <<- esc + 1; invokeRestart("muffleWarning") })
cat("escaped past expect_warning():", esc, "\n")
