# Review of lane mfxse, attack 1: J V J' on designs the lane did not use.
#   Rscript dev/review-mfxse-designs.R [lane|before] [design ...]
# Data seed 8801; Monte Carlo seed 8802, 2000 draws of the fit's joint normal.
#
# Three references per quantity, none of which reads get_coef(), set_coef()
# or get_vcov():
#   ref  the delta method over the WHOLE joint vector (every component the
#        joint covariance has: beta, betad, b, theta, ...), gradient by
#        central differences of get_predict() at estimates written into
#        fit$estimates by component; V = get_joint_cov(fit)$V (what fitted()
#        reads). A component the lane leaves out shows as mfx / ref < 1.
#   mc   SD over draws from N(est, V) of the same quantity.
#   fit  fitted(newdata)'s Est.Error (levels only).
a <- commandArgs(TRUE)
arm <- if (length(a)) a[1] else "lane"
only <- if (length(a) > 1) a[-1] else NULL
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(switch(arm, lane = c(lane, base, user), before = c(base, user),
                 c(arm, base, user)))
Sys.setenv(NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "arm:", arm, "\n")
quiet <- function(e) suppressWarnings(suppressMessages(e))
ns <- asNamespace("frmtmb")

set.seed(8801)
J <- 15; m <- 12; n <- J * m
d <- data.frame(g = factor(rep(seq_len(J), each = m)), x = rnorm(n),
                x2 = runif(n, -1, 1), t = rep(seq_len(m), J))
# h: a second grouping factor whose levels carry the same names as g's
d$h <- factor(sample(seq_len(J), n, TRUE), levels = seq_len(J))
d$a <- factor(sample(c("p", "q", "r"), n, TRUE))
d$bb <- factor(sample(c("u", "v"), n, TRUE))
d$g2 <- factor(sample(seq_len(J), n, TRUE), levels = seq_len(J))
d$w1 <- runif(n, 0.2, 0.8); d$w2 <- 1 - d$w1
u0 <- rnorm(J, 0, 0.7); u1 <- rnorm(J, 0, 0.4); uh <- rnorm(J, 0, 0.5)
uab <- rnorm(6, 0, 0.6)
ab <- as.integer(interaction(d$a, d$bb))
d$y <- 1 + 0.5 * d$x + u0[d$g] + u1[d$g] * d$x + uh[d$h] + uab[ab] +
  rnorm(n, 0, 0.6)
d$y2 <- -0.5 + 0.3 * d$x + 0.8 * u0[d$g] + rnorm(n, 0, 0.5)
d$ln <- exp(0.2 + 0.3 * d$x + u0[d$g] + rnorm(n, 0, exp(-0.5 + 0.5 * u1[d$g])))
d$expo <- runif(n, 0.5, 3)
d$cnt <- rpois(n, d$expo * exp(0.3 + 0.4 * d$x + u0[d$g]))
d$zc <- ifelse(runif(n) < plogis(-1 + uh[d$g]), 0L,
               rpois(n, exp(0.5 + 0.4 * d$x + u0[d$g])))
d$cat <- factor(apply(cbind(0, 0.3 + 0.5 * d$x + u0[d$g],
                            -0.2 - 0.4 * d$x + uh[d$g]) +
                        matrix(-log(-log(runif(3 * n))), n), 1, which.max))
d$wt <- runif(n, 0.5, 2)
d$ynl <- (2 + 0.4 * u0[d$g]) * exp(0.5 * d$x2) + rnorm(n, 0, 0.3)
d$xo <- d$x + rnorm(n, 0, 0.3); d$sx <- 0.3
d$ynull <- 1 + 0.5 * d$x + rnorm(n, 0, 0.6)   # no group effect: sd at 0

specs <- list(
  slope_cor = list(bf(y ~ x + (1 + x | g))),
  shared_lev = list(bf(y ~ x + (1 | g) + (1 | h))),
  inter_ab = list(bf(y ~ x + (1 | a:bb))),
  mm_w = list(bf(y ~ x + (1 | mm(g, g2, weights = cbind(w1, w2))))),
  sigma_g = list(bf(ln ~ x + (1 | g), sigma ~ 1 + (1 | g)), family = lognormal()),
  zi_g = list(bf(zc ~ x + (1 | g), zi ~ (1 | g)), family = zero_inflated_poisson()),
  categ = list(bf(cat ~ x + (1 | g)), family = categorical()),
  mv_id = list(mvbf(bf(y ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
                 set_rescor(FALSE)),
  offset = list(bf(cnt ~ x + offset(log(expo)) + (1 | g)), family = poisson()),
  weights = list(bf(y | weights(wt) ~ x + (1 | g))),
  nlpar = list(bf(ynl ~ alpha * exp(kk * x2), alpha ~ 1 + (1 | g), kk ~ 1, nl = TRUE)),
  me = list(bf(y ~ me(xo, sx) + (1 | g))),
  latent_ar = list(bf(cnt ~ x + (1 | g) + ar(time = t, gr = g, cov = TRUE)), family = poisson()),
  reml = list(bf(y ~ x + (1 + x | g)), REML = TRUE),
  singular = list(bf(ynull ~ x + (1 | g)))
)
if (!is.null(only)) specs <- specs[only]

# the joint vector, component by component, as estimates hold it
jv_map <- function(fit, rn) {
  est <- fit$estimates
  out <- list()
  for (cp in unique(rn)) {
    k <- sum(rn == cp)
    e <- est[[cp]]
    pos <- seq_along(e)
    if (cp == "betad" && length(fit$frame[["betad_fixed_idx"]])) {
      pos <- setdiff(pos, fit$frame[["betad_fixed_idx"]])
    }
    if (length(pos) != k) return(paste("unmapped", cp, k, length(e)))
    out[[cp]] <- pos
  }
  out
}
jv_get <- function(fit, rn, mp) {
  v <- numeric(length(rn))
  for (cp in names(mp)) v[rn == cp] <- fit$estimates[[cp]][mp[[cp]]]
  v
}
jv_put <- function(fit, rn, mp, v) {
  for (cp in names(mp)) fit$estimates[[cp]][mp[[cp]]] <- v[rn == cp]
  fit
}

rows_of <- function(fit, nd, ...) {
  marginaleffects::get_predict(fit, newdata = nd, ...)$estimate
}

run_design <- function(k, s) {
  args <- c(list(s[[1L]], data = d), s[-1L])
  fit <- tryCatch(quiet(do.call(frm, args)), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-11s FITERR %s\n", k, conditionMessage(fit)))
    return(invisible())
  }
  conv <- fit$opt$convergence
  lost <- names(ns$sdr_of(fit)$se_lost)
  jc <- ns$get_joint_cov(fit)
  rn <- jc$names; V <- as.matrix(jc$V)
  mp <- jv_map(fit, rn)
  if (is.character(mp)) { cat(k, mp, "\n"); return(invisible()) }
  v0 <- jv_get(fit, rn, mp)
  resp <- if (k == "mv_id") "y" else NULL
  nd <- d[c(5, 61, 122, 170), ]
  # the quantities, each a function of a fit
  qs <- list(
    pred = function(f) rows_of(f, nd, resp = resp),
    avg = function(f) {
      p <- marginaleffects::get_predict(f, newdata = d, resp = resp)
      if (!is.null(p$group)) tapply(p$estimate, p$group, mean) else mean(p$estimate)
    },
    slope = function(f) {
      h <- 1e-4
      hi <- d; hi$x <- hi$x + h / 2; lo <- d; lo$x <- lo$x - h / 2
      ph <- marginaleffects::get_predict(f, newdata = hi, resp = resp)
      pl <- marginaleffects::get_predict(f, newdata = lo, resp = resp)
      dd <- (ph$estimate - pl$estimate) / h
      if (!is.null(ph$group)) tapply(dd, ph$group, mean) else mean(dd)
    })
  if (k == "nlpar") qs$slope <- NULL
  ws <- character(0)
  catch <- function(e) withCallingHandlers(e, warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 70))
    invokeRestart("muffleWarning")
  })
  mfx <- list(
    pred = catch(tryCatch(predictions(fit, newdata = nd, resp = resp),
                          error = function(e) conditionMessage(e))),
    avg = catch(tryCatch(avg_predictions(fit, resp = resp),
                         error = function(e) conditionMessage(e))),
    slope = if (!is.null(qs$slope)) catch(tryCatch(
      avg_slopes(fit, variables = "x", resp = resp),
      error = function(e) conditionMessage(e))))
  ftse <- tryCatch(quiet(fitted(fit, newdata = nd, resp = resp))[, "Est.Error"],
                   error = function(e) NA)
  # Monte Carlo draws of the joint normal
  set.seed(8802)
  ev <- eigen((V + t(V)) / 2, symmetric = TRUE)
  R <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(ev$values))
  Z <- matrix(rnorm(2000 * length(v0)), ncol = length(v0))
  draws <- sweep(Z %*% t(R), 2, v0, "+")
  ncp <- length(get_coef(fit))
  cat(sprintf("\n== %s  conv %d  joint rows %d  get_coef %d  lost: %s  min eig V %.2e\n",
              k, conv, length(rn), ncp, paste(lost, collapse = ","),
              min(ev$values)))
  for (q in names(qs)) {
    f0 <- qs[[q]](fit)
    G <- matrix(0, length(f0), length(v0))
    for (j in seq_along(v0)) {
      hj <- 1e-5 * max(1, abs(v0[j]))
      vp <- v0; vp[j] <- vp[j] + hj
      vm <- v0; vm[j] <- vm[j] - hj
      G[, j] <- (qs[[q]](jv_put(fit, rn, mp, vp)) -
                   qs[[q]](jv_put(fit, rn, mp, vm))) / (2 * hj)
    }
    ref <- sqrt(pmax(rowSums((G %*% V) * G), 0))
    nmc <- if (q == "pred") 2000 else 400
    sims <- t(vapply(seq_len(nmc), function(i) {
      as.numeric(qs[[q]](jv_put(fit, rn, mp, draws[i, ])))
    }, numeric(length(f0))))
    if (length(f0) == 1) sims <- t(sims)
    mc <- apply(sims, 2, stats::sd)
    mo <- mfx[[q]]
    if (is.character(mo)) {
      cat(sprintf("  %-5s ERROR %s\n", q, mo)); next
    }
    se <- mo$std.error
    est_d <- max(abs(mo$estimate - as.numeric(f0)))
    cat(sprintf("  %-5s mfx/ref %s  | mfx/mc %s  | est diff %.1e\n", q,
                paste(sprintf("%.4f", se / ref), collapse = " "),
                paste(sprintf("%.3f", se / mc), collapse = " "), est_d))
    if (q == "pred") {
      cat(sprintf("        mfx/fitted %s  | ref/fitted %s\n",
                  paste(sprintf("%.4f", se / ftse), collapse = " "),
                  paste(sprintf("%.4f", ref / ftse), collapse = " ")))
    }
  }
  if (length(ws)) cat("  warnings:", paste(unique(ws), collapse = " || "), "\n")
}
for (k in names(specs)) {
  r <- tryCatch(run_design(k, specs[[k]]), error = function(e) {
    cat(sprintf("%-11s RUNERR %s\n", k, conditionMessage(e)))
  })
}
cat("done\n")
