# Review of lane parity2, item 5300: variables() of hand-built draws of a
# latent ar model, arm base|lane, and how err[k] enters the shift of the
# linear predictor under cov = FALSE (regression on err[n], err[prev(n)]).
#   Rscript dev/reviews/parity2-review-err2.R base|lane
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-parity2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
q <- function(e) suppressMessages(suppressWarnings(e))
set.seed(4)
dp <- data.frame(t = rep(1:8, 6), g = factor(rep(1:6, each = 8)),
                 x = rnorm(48))
e <- as.vector(replicate(6, arima.sim(list(ar = 0.5), 8, sd = 0.6)))
dp$y <- rpois(48, exp(0.5 + 0.3 * dp$x + e))
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
fit <- q(frm(bf(y ~ x + ar(time = t, gr = g)), family = poisson(), data = dp))
ds <- hand_draws(fit)
v <- variables(ds)
cat("variables (no err/zerr):",
    paste(v[!grepl("^(err|zerr)\\[", v)], collapse = " "), "\n")
cat("count", length(v), " zerr", sum(grepl("^zerr\\[", v)), " err",
    sum(grepl("^err\\[", v)), "\n")
if (arm == "lane") {
  set.seed(5)
  sh <- dp[sample(48), ]
  rownames(sh) <- NULL
  fs <- q(frm(bf(y ~ x + ar(time = t, gr = g)), family = poisson(), data = sh))
  lp <- fitted(fs, scale = "linear")[, "Estimate"]
  fe <- fixef(fs)[, "Estimate"]
  s <- lp - (fe[["Intercept"]] + fe[["x"]] * sh$x)
  bk <- Filter(function(b) identical(b$covstruct, "latres"),
               fs$frame$re_blocks)[[1]]
  err <- frmtmb::latres_err(fs$estimates$theta[bk$theta_idx],
                            fs$estimates$b[bk$b_idx], bk)
  prev <- vapply(seq_len(nrow(sh)), function(i) {
    j <- which(sh$g == sh$g[i] & sh$t == sh$t[i] - 1)
    if (length(j)) j else NA_integer_
  }, 1L)
  lag1 <- ifelse(is.na(prev), 0, err[pmax(prev, 1L)])
  m <- lm(s ~ 0 + err + lag1)
  cat("shuffled rows, shift ~ err[n] + err[prev(n)]: coef",
      format(coef(m), digits = 6), " max |resid|", max(abs(resid(m))), "\n")
  o <- order(sh$g, sh$t)
  e2 <- err[order(o)]
  lag2 <- ifelse(is.na(prev), 0, e2[pmax(prev, 1L)])
  m2 <- lm(s ~ 0 + e2 + lag2)
  cat("same with err permuted as (g, t) order: max |resid|",
      max(abs(resid(m2))), "\n")
}
