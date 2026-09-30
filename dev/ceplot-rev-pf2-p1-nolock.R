# Copy of dev/postfit2-rev-p1-nolock.R for the ceplot review, with the library
# paths of this round (wt-ceplot-lib, rellib-r4); otherwise unchanged.
# Reviewer, punch round 1: what ce_locked_vars() protects. P1 of
# postfit2-rev-p1-attack.R (y ~ x + trt + (1 | trt:subj), data seed 41),
# trt = "b" set, subj unset, before and after replacing
# ce_locked_vars() by a function that locks nothing (in process).
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib", "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
set.seed(41)
d1 <- expand.grid(rep = 1:4, trt = factor(c("a", "b")), subj = factor(1:12))
d1$x <- rnorm(nrow(d1))
d1$y <- rnorm(nrow(d1), 1 + 0.5 * d1$x + (d1$trt == "b") +
                rnorm(24, 0, 1)[interaction(d1$trt, d1$subj)], 0.5)
f1 <- frm(bf(y ~ x + trt + (1 | trt:subj)), family = gaussian(), data = d1)
tpl <- f1$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) f1$estimates[[cp]]))
set.seed(1)
M <- matrix(rep(est, each = 200) + rnorm(200 * length(est), 0, 0.05), 200,
            dimnames = list(NULL, frmtmb::brms_par_labels(f1)))
ds <- structure(list(stanfit = NULL, draws = cbind(frmtmb.sample:::draws_to_natural(M, f1), lp__ = 0),
                     fit = f1), class = "frmtmb_draws")
go <- function(lab) {
  w <- suppressMessages(conditional_effects(f1, "x", resolution = 3, re_formula = NULL, conditions = list(trt = "b")))$x
  b <- suppressMessages(conditional_effects(f1, "x", resolution = 3, re_formula = NULL, conditions = list(trt = "b"), band = "boot", boot = 40, seed = 5))$x
  d <- suppressMessages(conditional_effects(ds, "x", resolution = 3, re_formula = NULL, conditions = list(trt = "b"), seed = 1))$x
  cat(lab, "wald est", f3(w$estimate__), "| boot mean of refits", f3(colMeans(attr(suppressMessages(conditional_effects(f1, "x", resolution = 3, re_formula = NULL, conditions = list(trt = "b"), band = "boot", boot = 40, seed = 5)), "boot")$t)),
      "| draws est", f3(d$estimate__), "\n")
}
go("lane:          ")
assignInNamespace("ce_locked_vars", function(fit) character(0), "frmtmb")
go("mutant nolock: ")
