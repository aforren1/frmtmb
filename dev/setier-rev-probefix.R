# Reviewer of lane setier: a candidate fix for se_curvature_real(),
# patched into the lane's namespace for this process only. The probe
# predicts a loss of 2 * grad_tol each way; the fix also takes half the
# step and keeps a direction only when the loss scales as a quadratic
# (ratio 3 to 5.5 for twice the step), so a straight step off a curved
# ridge, which loses as the fourth power, is not taken as curvature.
#   Rscript dev/setier-rev-probefix.R
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
fixed <- function(fit, p, free, dir_free, lambda) {
  tol <- fit$control$grad_tol %||% 1e-3
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  d <- numeric(length(p))
  d[free] <- dir_free * sqrt(4 * tol / lambda)
  f0 <- tryCatch(obj$fn(p), error = function(e) NA_real_)
  fu <- tryCatch(obj$fn(p + d), error = function(e) NA_real_)
  fd <- tryCatch(obj$fn(p - d), error = function(e) NA_real_)
  fu2 <- tryCatch(obj$fn(p + d / 2), error = function(e) NA_real_)
  fd2 <- tryCatch(obj$fn(p - d / 2), error = function(e) NA_real_)
  l <- c(fu, fd) - f0
  r <- l / (c(fu2, fd2) - f0)
  # a quadratic loses four times as much at twice the step; a straight
  # step off a curved ridge loses as the fourth power (ratio 16)
  all(is.finite(c(f0, fu, fd, fu2, fd2))) && all(l >= tol) &&
    all(r > 3 & r < 5.5)
}
environment(fixed) <- ns
probe_log <- character()
for (arm in c("lane", "fixed")) {
  if (arm == "fixed") {
    unlockBinding("se_curvature_real", ns)
    assign("se_curvature_real", fixed, envir = ns)
  }
  cat("\n== arm", arm, "\n")
  # c0k, seed 77
  set.seed(77)
  G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(G, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  f <- suppressMessages(suppressWarnings(frm(bf(yn ~ c0 + exp(a)^k,
    c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE), data = dn)))
  lost <- ns$sdr_of(f)$se_lost
  cat("c0k seed 77: lost", paste(names(lost), collapse = ","), "\n")
  # raw polynomials, degree 5 and 6, with and without (1 | g)
  for (deg in 5:6) for (re in c(FALSE, TRUE)) {
    set.seed(deg)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
    fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                           if (re) "+ (1 | g)"))
    f <- suppressMessages(suppressWarnings(frm(fo, data = d)))
    se <- suppressWarnings(fixef(f))[, "Est.Error"]
    ref <- if (re) {
      sqrt(diag(as.matrix(vcov(suppressMessages(suppressWarnings(
        lmer(fo, data = d, REML = FALSE)))))))
    } else sqrt(diag(vcov(lm(fo, data = d))) * (200 - deg - 1) / 200)
    cat(sprintf("poly degree %d re %-5s: lost %s | max |SE/ref - 1| %.3g\n",
                deg, re, length(ns$sdr_of(f)$se_lost),
                max(abs(unname(se) / unname(ref) - 1))))
  }
}
