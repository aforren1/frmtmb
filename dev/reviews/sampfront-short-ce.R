# Short check: what conditional_effects(method = "posterior_predict",
# re_formula = NULL) gives on laplace draws of bf(y ~ x, sigma ~ (1 | h)),
# with its warnings shown.   Rscript dev/reviews/sampfront-short-ce.R [base]
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", if (base) "base" else "lane", "\n")
lap_pair <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  lap <- structure(list(stanfit = NULL, draws = M, fit = fit),
                   class = "frmtmb_draws")
  lap$draws <- M[, setdiff(colnames(M), inner), drop = FALSE]
  lap
}
set.seed(5)
n <- 60
d <- data.frame(x = rnorm(n), h = factor(rep(1:5, each = 12)))
d$y <- 0.3 * d$x + rnorm(n, 0, exp(rnorm(5, 0, 0.5))[d$h])
fit <- suppressWarnings(frm(bf(y ~ x, sigma ~ (1 | h)), data = d))
lap <- lap_pair(fit)
ws <- character()
r <- withCallingHandlers(
  tryCatch(conditional_effects(lap, effects = "x", re_formula = NULL,
             method = "posterior_predict",
             conditions = data.frame(h = factor(2, levels = 1:5))),
           error = function(e) e),
  warning = function(w) { ws <<- c(ws, conditionMessage(w))
                          invokeRestart("muffleWarning") })
if (inherits(r, "error")) cat("REFUSED:", conditionMessage(r), "\n") else {
  print(head(r[[1]][, c("x", "estimate__", "lower__", "upper__")], 3))
}
cat("warnings:", length(ws), "\n"); print(table(ws))
