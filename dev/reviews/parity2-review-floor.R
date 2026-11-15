# Review of lane parity2: the lane's frmtmb.sample on the BASE core
# (rellib-r21, 0.82.0), which lacks ce_snl() and latres_err(): does it
# fail loudly where it needs them, and answer where it does not?
#   Rscript dev/reviews/parity2-review-floor.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
loadNamespace("frmtmb", lib.loc = "C:/Users/adf44/source/r/rellib-r21")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample, lib.loc = "C:/Users/adf44/source/r/wt-parity2-lib")
})
cat("frmtmb", find.package("frmtmb"), " frmtmb.sample",
    find.package("frmtmb.sample"), "\n")
cat("frmtmb.sample DESCRIPTION Imports frmtmb:",
    grep("frmtmb \\(", packageDescription("frmtmb.sample")$Imports,
         value = TRUE), "\n")
hand_draws <- function(fit, n = 6L, seed = 1L, sd = 0.02) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, sd),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit, data = fit$data),
            class = "frmtmb_draws")
}
tr <- function(label, e) {
  r <- tryCatch(e, error = function(err) err)
  cat(label, ": ", if (inherits(r, "error")) paste("ERROR", conditionMessage(r))
      else "answers", "\n", sep = "")
}
set.seed(11)
dg <- data.frame(x = rnorm(48), g = factor(rep(1:6, each = 8)))
dg$y <- 1 + 0.5 * dg$x + rnorm(6, sd = 0.8)[dg$g] + rnorm(48, sd = 0.5)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dg)
ds <- hand_draws(fit)
tr("ce re_formula = NA", conditional_effects(ds, "x"))
tr("ce re_formula = NULL", conditional_effects(ds, "x", re_formula = NULL))
tr("loo_epred subset", suppressWarnings(suppressMessages(
  loo_epred(ds, newdata = dg[1:5, ]))))
set.seed(4)
dp <- data.frame(t = rep(1:8, 6), g = factor(rep(1:6, each = 8)),
                 x = rnorm(48))
dp$y <- rpois(48, exp(0.5 + 0.3 * dp$x))
fa <- suppressWarnings(suppressMessages(
  frm(bf(y ~ x + ar(time = t, gr = g)), family = poisson(), data = dp)))
tr("latres draws", hand_draws(fa))
