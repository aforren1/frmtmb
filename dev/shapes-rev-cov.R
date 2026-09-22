# Reviewer, priority 2: an INDEPENDENT coverage measurement of
# predict()'s interval, built without reading dev/shapes-coverage.R's
# arms, plus the two things that script has no arm for.
#
# What is different here on purpose:
#
#  * a POSITIVE control. The lane's harness carries only a control that
#    must FAIL. On the gaussian design the exact frequentist prediction
#    interval is known in closed form,
#      xhat'b +- t(n - p, 0.975) * s * sqrt(1 + x'(X'X)^-1 x),
#    and its coverage is 0.95 by construction, so an arm that reports
#    anything else says the harness is wrong rather than predict().
#  * the mixed design's ATTRIBUTION is tested, not asserted. The lane
#    says the shortfall is the omitted Var(b | y). This adds that term
#    to the draws and reports the coverage again. The conditional
#    variance of a random intercept is
#      sigma^2 tau^2 / (sigma^2 + n_g tau^2),
#    computed here from the fit's own estimates rather than from any
#    frmtmb accessor, so it is an outside number.
#  * the seed streams of the arms are DISJOINT by construction and the
#    script prints the number of distinct seeds it used, so a reused
#    stream would show as a count.
#
#   Rscript dev/shapes-rev-cov.R <design> <nrep>
# designs: gauss, mixed

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

a <- commandArgs(trailingOnly = TRUE)
design <- if (length(a)) a[1] else "gauss"
nrep <- if (length(a) > 1) as.integer(a[2]) else 220L
n <- 60L; m <- 10L; ndraws <- 1000L
SEED0 <- 771000L          # deliberately NOT the lane's 20260917

seeds_used <- integer(0)
sset <- function(s) { seeds_used <<- c(seeds_used, s); set.seed(s) }

hits <- list(); est_err <- list()
add <- function(nm, v) hits[[nm]] <<- c(hits[[nm]], v)
fail <- 0L; tot <- 0L

for (r in seq_len(nrep)) {
  s <- SEED0 + 100L * r          # 100 apart: no arm's stream can meet
  sset(s)                         # another replicate's
  if (design == "gauss") {
    d <- data.frame(x = rnorm(n))
    d$y <- rnorm(n, -0.4 + 1.3 * d$x, 0.8)
    form <- bf(y ~ x) + gaussian()
    nd <- data.frame(x = rnorm(m))
    ynew <- rnorm(m, -0.4 + 1.3 * nd$x, 0.8)
  } else {
    ng <- 12L
    d <- data.frame(x = rnorm(n),
                    g = factor(rep(seq_len(ng), length.out = n)))
    u <- rnorm(ng, 0, 0.7)
    d$y <- rnorm(n, -0.4 + 1.3 * d$x + u[d$g], 0.8)
    form <- bf(y ~ x + (1 | g)) + gaussian()
    gsel <- sample(seq_len(ng), m, TRUE)
    nd <- data.frame(x = rnorm(m),
                     g = factor(gsel, levels = levels(d$g)))
    ynew <- rnorm(m, -0.4 + 1.3 * nd$x + u[gsel], 0.8)
  }
  fit <- tryCatch(suppressWarnings(frm(form, data = d)),
                  error = function(e) NULL)
  if (is.null(fit)) { fail <- fail + 1L; next }

  sset(s + 1L)
  dj <- tryCatch(suppressWarnings(
    predict(fit, newdata = nd, ndraws = ndraws, summary = FALSE)),
    error = function(e) NULL)
  sset(s + 2L)
  dp <- tryCatch(suppressWarnings(
    predict(fit, newdata = nd, ndraws = ndraws, summary = FALSE,
            param_uncertainty = FALSE)), error = function(e) NULL)
  if (is.null(dj) || is.null(dp)) { fail <- fail + 1L; next }
  qi <- function(D) t(apply(D, 2L, stats::quantile, c(0.025, 0.975),
                            names = FALSE))
  Qj <- qi(dj); Qp <- qi(dp)
  add("joint",  sum(ynew >= Qj[, 1] & ynew <= Qj[, 2]))
  add("plugin", sum(ynew >= Qp[, 1] & ynew <= Qp[, 2]))

  # the POSITIVE control: the textbook prediction interval, computed
  # from lm() so that no frmtmb code is on its path at all
  if (design == "gauss") {
    lmf <- stats::lm(y ~ x, data = d)
    pi <- stats::predict(lmf, newdata = nd, interval = "prediction",
                         level = 0.95)
    add("exact", sum(ynew >= pi[, "lwr"] & ynew <= pi[, "upr"]))
  } else {
    # lme4's own predictive interval has no closed form; the positive
    # control on this design is the gaussian one above
    add("exact", NA_integer_)
  }

  # the attribution test: add the conditional variance of the random
  # intercept to the joint draws and count again
  if (design == "mixed") {
    sg <- exp(as.numeric(fit$estimates[["betad"]][1L]))   # log sigma
    th <- as.numeric(fit$estimates[["theta"]][1L])        # log tau
    tau <- exp(th)
    ng_obs <- as.integer(table(d$g))
    cv <- sg^2 * tau^2 / (sg^2 + ng_obs * tau^2)
    sd_add <- sqrt(cv[as.integer(nd$g)])
    sset(s + 3L)
    dj2 <- dj + matrix(stats::rnorm(length(dj), 0,
                                    rep(sd_add, each = nrow(dj))),
                       nrow(dj), ncol(dj))
    Q2 <- qi(dj2)
    add("joint_plus_condvar",
        sum(ynew >= Q2[, 1] & ynew <= Q2[, 2]))
    est_err[[length(est_err) + 1L]] <-
      c(sigma = sg, tau = tau, mean_condsd = mean(sd_add))
  }

  # is the interval for a NEW OBSERVATION? On a gaussian fit the
  # identity is Var(pred) = Var(fitted) + sigma^2 (+ the parameter
  # draw's own share). Report the ratio; it is a check on the arm, not
  # a coverage count.
  if (design == "gauss" && r <= 20L) {
    fv <- fitted(fit, newdata = nd)
    sg <- exp(as.numeric(fit$estimates[["betad"]][1L]))
    got <- apply(dp, 2L, stats::sd)
    want <- sqrt(sg^2)            # plug-in: no parameter draw at all
    est_err[[length(est_err) + 1L]] <-
      c(ratio_plugin_sd_over_sigma = max(abs(got / want - 1)),
        fitted_se_max = max(fv[, "Est.Error"]))
  }
  tot <- tot + m
}

cat("design:", design, " replicates:", length(hits$joint),
    " failed fits:", fail, "\n")
cat("n:", n, " new points:", m, " ndraws:", ndraws, " SEED0:", SEED0,
    "\n")
cat("distinct seeds used:", length(unique(seeds_used)), "of",
    length(seeds_used), "set.seed calls\n")
cat("script: dev/shapes-rev-cov.R\n\n")
for (nm in names(hits)) {
  v <- hits[[nm]]
  if (all(is.na(v))) { cat(sprintf("%-20s n/a\n", nm)); next }
  k <- sum(v)
  ci <- stats::binom.test(k, tot, 0.95)$conf.int
  pr <- v / m
  cat(sprintf("%-20s %5d of %5d = %.4f  (%.4f, %.4f) binomial; ",
              nm, k, tot, k / tot, ci[1], ci[2]))
  cat(sprintf("replicate mean %.4f se %.4f  z=%.2f\n", mean(pr),
              stats::sd(pr) / sqrt(length(pr)),
              (mean(pr) - 0.95) / (stats::sd(pr) / sqrt(length(pr)))))
}
if (length(est_err)) {
  E <- do.call(rbind, est_err)
  cat("\nside measurements (column means over ", nrow(E), " rows)\n",
      sep = "")
  print(round(colMeans(E), 5))
}
