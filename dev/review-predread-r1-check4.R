# Short check, item 4: predict(newdata) and posterior_predict(newdata) under cov = TRUE
# against brms at the same values, and log_lik(newdata =). Seed 8201. Also run on the
# merged tree (REV_LIB=<merge lib>) for a latent-residual family.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(8201)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(rep(sprintf("g%02d", 1:30), each = 10)))
d$t <- rep(1:10, 30)
e <- unlist(lapply(1:30, function(i) as.numeric(arima.sim(list(ar = 0.6), 10, sd = 0.7))))
d$y <- 0.5 + 0.8 * d$x + e
d$yt <- 0.5 + 0.8 * d$x + e * 1.1
le <- unlist(lapply(1:30, function(i) as.numeric(arima.sim(list(ar = 0.6), 10, sd = 0.5))))
d$cnt <- rpois(n, exp(0.5 + 0.4 * d$x + le))
d$bin <- rbinom(n, 1, plogis(0.3 + 0.7 * d$x + le))
nd <- data.frame(x = c(-1, 0, 1, 0.5), t = c(1L, 2L, 3L, 4L), g = factor("g01", levels = levels(d$g)))
fmt <- function(v) paste(sprintf("%.3f", v), collapse = " ")
msg <- function(e) tryCatch(quiet(e), error = function(err) structure(NA, err = substr(conditionMessage(err), 1, 170)))
specs <- list(
  ar_cov = c("bf(y ~ x + ar(time = t, gr = g, cov = TRUE))", "gaussian()"),
  ma_cov = c("bf(y ~ x + ma(time = t, gr = g, cov = TRUE))", "gaussian()"),
  arma_cov = c("bf(y ~ x + arma(time = t, gr = g, cov = TRUE))", "gaussian()"),
  cosy = c("bf(y ~ x + cosy(time = t, gr = g))", "gaussian()"),
  student_ar = c("bf(y ~ x + ar(time = t, gr = g, cov = TRUE))", "student()"),
  sigma_ar = c("bf(y ~ x + ar(time = t, gr = g, cov = TRUE), sigma ~ x)", "gaussian()"),
  poisson_ar = c("bf(cnt ~ x + ar(time = t, gr = g, cov = TRUE))", "poisson()"),
  bernoulli_ar = c("bf(bin ~ x + ar(time = t, gr = g, cov = TRUE))", "bernoulli()"))
for (nm in names(specs)) {
  cat("\n===", nm, specs[[nm]], "\n")
  pr <- tryCatch(rev_pair(specs[[nm]][1], specs[[nm]][2], d, ndraws = 4000), error = function(e) e)
  fit <- NULL
  fx <- NULL
  if (inherits(pr, "error")) {
    cat("pair ERR:", substr(conditionMessage(pr), 1, 170), "\n")
    fit <- msg(frm(rev_ev(specs[[nm]][1], "frmtmb"), family = rev_ev(specs[[nm]][2], "frmtmb"), data = d))
    if (!is.null(attr(fit, "err"))) {
      cat("frm ERR:", attr(fit, "err"), "\n")
      next
    }
  } else {
    fit <- pr$fit
    fx <- pr$fx
  }
  if (!is.null(fx)) {
    bp <- posterior_predict(fx, newdata = nd)
    cat("brms pp      mean", fmt(colMeans(bp)), " sd", fmt(apply(bp, 2, sd)), " cor", fmt(cor(bp)[1, 2:4]), "\n")
    be <- posterior_epred(fx, newdata = nd)
    cat("brms epred   mean", fmt(colMeans(be)), " sd", fmt(apply(be, 2, sd)), "\n")
  }
  set.seed(5)
  fp <- msg(predict(fit, newdata = nd, summary = FALSE, ndraws = 4000, propagate_error = FALSE))
  if (is.null(attr(fp, "err"))) {
    cat("predict      mean", fmt(colMeans(fp)), " sd", fmt(apply(fp, 2, sd)), " cor", fmt(cor(fp)[1, 2:4]), "\n")
  } else cat("predict ERR", attr(fp, "err"), "\n")
  fv <- msg(fitted(fit, newdata = nd))
  if (is.null(attr(fv, "err"))) cat("fitted       ", fmt(fv[, "Estimate"]), " se", fmt(fv[, "Est.Error"]), "\n") else cat("fitted ERR", attr(fv, "err"), "\n")
  fm <- msg(t(vapply(simulate(fit, nsim = 4000, seed = 3, newdata = nd), as.numeric, numeric(4))))
  if (is.null(attr(fm, "err"))) cat("simulate     mean", fmt(colMeans(fm)), " sd", fmt(apply(fm, 2, sd)), " cor", fmt(cor(fm)[1, 2:4]), "\n") else cat("simulate ERR", attr(fm, "err"), "\n")
  ds <- rev_draws_at(fit, 4000L)
  set.seed(6)
  dp <- msg(frmtmb.sample::posterior_predict(ds, newdata = nd))
  if (is.null(attr(dp, "err"))) cat("draws pp     mean", fmt(colMeans(dp)), " sd", fmt(apply(dp, 2, sd)), " cor", fmt(cor(dp)[1, 2:4]), "\n") else cat("draws pp ERR", attr(dp, "err"), "\n")
  de <- msg(frmtmb.sample::posterior_epred(ds, newdata = nd))
  if (is.null(attr(de, "err"))) cat("draws epred  mean", fmt(colMeans(de)), " sd", fmt(apply(de, 2, sd)), "\n") else cat("draws epred ERR", attr(de, "err"), "\n")
  yn <- all.vars(rev_ev(specs[[nm]][1], "frmtmb")$formula)[1L]
  nl <- nd
  nl[[yn]] <- if (yn == "bin") c(0L, 1L, 1L, 0L) else if (yn == "cnt") c(1L, 2L, 3L, 1L) else c(0.1, 0.2, 1, 0.5)
  if (!is.null(fx)) {
    bl <- msg(log_lik(fx, newdata = nl))
    cat("brms log_lik newdata:", if (is.null(attr(bl, "err"))) fmt(colMeans(bl)) else attr(bl, "err"), "\n")
  }
  fl <- msg(frmtmb.sample::log_lik(ds, newdata = nl))
  cat("draws log_lik newdata:", if (is.null(attr(fl, "err"))) fmt(colMeans(fl)) else attr(fl, "err"), "\n")
}
