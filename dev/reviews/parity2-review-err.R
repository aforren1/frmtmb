# Review of lane parity2, item 5300: latres_err() against the fit's own
# shift of the linear predictor, per DATA row, on rows that are not in
# (g, t) order; and the err names of frmtmb.sample on ma, cov = TRUE and
# a multivariate latent model (hand-built draws).
#   Rscript dev/reviews/parity2-review-err.R
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
q <- function(e) suppressMessages(suppressWarnings(e))
set.seed(4)
dp <- data.frame(t = rep(1:8, 6), g = factor(rep(1:6, each = 8)),
                 x = rnorm(48))
e <- as.vector(replicate(6, arima.sim(list(ar = 0.5), 8, sd = 0.6)))
dp$y <- rpois(48, exp(0.5 + 0.3 * dp$x + e))
dp$y2 <- rpois(48, exp(0.2 + 0.1 * dp$x))
set.seed(5)
sh <- dp[sample(48), ]
rownames(sh) <- NULL
shift_of <- function(fit, data) {
  lp <- fitted(fit, scale = "linear")[, "Estimate"]
  fe <- fixef(fit)[, "Estimate"]
  xb <- fe[["Intercept"]] + fe[["x"]] * data$x
  lp - xb
}
blk_err <- function(fit) {
  bk <- Filter(function(b) identical(b$covstruct, "latres"),
               fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  b <- fit$estimates$b[bk$b_idx]
  list(err = frmtmb::latres_err(th, b, bk), bk = bk, th = th)
}
for (cv in c(TRUE, FALSE)) {
  fit <- q(frm(bf(y ~ x + ar(time = t, gr = g, cov = cv)), family = poisson(),
               data = sh))
  s <- shift_of(fit, sh)
  be <- blk_err(fit)
  if (cv) {
    cat("cov = TRUE, shuffled rows: max |err - shift| =",
        max(abs(be$err - s)), " max |shift| =", max(abs(s)), "\n")
    o <- order(sh$g, sh$t)
    cat("   against shift in (g, t) order:", max(abs(be$err - s[o])), "\n")
  } else {
    arv <- fit$estimates$theta[be$bk$theta_idx]
    cat("cov = FALSE theta of the block:", format(arv, digits = 4), "\n")
    prev <- vapply(seq_len(nrow(sh)), function(i) {
      j <- which(sh$g == sh$g[i] & sh$t == sh$t[i] - 1)
      if (length(j)) j else NA_integer_
    }, 1L)
    # with one ar coefficient, the shift of row n is ar * err[prev(n)]
    ok <- !is.na(prev)
    ratio <- s[ok] / be$err[prev[ok]]
    cat("cov = FALSE shuffled: shift / err[prev] range",
        format(range(ratio), digits = 6), " shift at first times",
        format(max(abs(s[!ok])), digits = 3), "\n")
    o <- order(sh$g, sh$t)
    e2 <- be$err[o]
    ratio2 <- s[ok] / e2[prev[ok]]
    cat("   if err were in (g, t) order: range", format(range(ratio2),
                                                         digits = 6), "\n")
  }
}

hand_draws <- function(fit, n = 6L, seed = 1L, sd = 0.02) {
  tpl <- fit$frame[["par_template"]]
  fx <- fit$frame[["betad_fixed_idx"]]
  est <- unlist(lapply(names(tpl), function(cp) {
    v <- fit$estimates[[cp]]
    if (cp == "betad" && length(fx)) v <- v[-fx]
    v
  }))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, sd),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit, data = fit$data),
            class = "frmtmb_draws")
}
show <- function(label, fit) {
  ds <- tryCatch(hand_draws(fit), error = function(e) e)
  if (inherits(ds, "error")) {
    cat(label, ": ERROR ", conditionMessage(ds), "\n", sep = "")
    return(invisible())
  }
  v <- variables(ds)
  cat(label, ": n", length(v), " first err at", grep("^err", v)[1],
      " non-err-tail: ", paste(v[!grepl("^err", v) &
                                    !grepl("^(zerr|latent_state)\\[", v)],
                               collapse = " "), "\n", sep = "")
  cat("   err names:", paste(head(grep("^err", v, value = TRUE), 2),
                             collapse = " "), "... count",
      length(grep("^err", v)), "\n")
}
show("ar cov FALSE", q(frm(bf(y ~ x + ar(time = t, gr = g)),
                          family = poisson(), data = dp)))
show("ar cov TRUE", q(frm(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)),
                          family = poisson(), data = dp)))
show("ma", q(frm(bf(y ~ x + ma(time = t, gr = g)), family = poisson(),
                 data = dp)))
show("mv latent", q(frm(bf(y ~ x + ar(time = t, gr = g)) + bf(y2 ~ x) +
                          set_rescor(FALSE), family = poisson(), data = dp)))
