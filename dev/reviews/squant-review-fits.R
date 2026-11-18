# Review of lane squant, item 3: independent qgam comparisons (two
# parametric, one smooth), rq() at small smoothing, and the bias claim
# "mu estimates the smoothed quantile, bias ~ -(pi^2/6)(lambda sigma)^2
# f'(q)/f(q)" on skewed residuals (population by integration, and one
# large-sample fit).
#
# Run (worktree root):
#   Rscript dev/reviews/squant-review-fits.R > dev/reviews/squant-review-fits.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(qgam)
  library(quantreg)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "qgam", format(packageVersion("qgam")), "\n")
quiet <- function(expr) {
  ws <- character()
  v <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 100))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = ws)
}
same_model <- function(qf) {
  co <- unique(qf$family$getCo())
  th <- qf$family$getTheta()
  list(sigma = exp(th), smoothing = co / exp(th), co = co)
}
qgam_ll <- function(qf, y) {
  -qf$family$aic(y, stats::fitted(qf), qf$family$getTheta(),
                 rep(1, length(y)), NULL) / 2
}

## 3a. parametric cases, not the lane's data ------------------------------
cat("\n== 3a. qgam, parametric, sigma held at exp(lsig)\n")
set.seed(301)
n <- 500
d1 <- data.frame(x1 = runif(n), x2 = rnorm(n), f = factor(sample(letters[1:3],
                                                                  n, TRUE)))
d1$y <- 1 + 2 * d1$x1 - 0.5 * d1$x2 + c(a = 0, b = 0.7, c = -0.4)[d1$f] +
  (0.4 + d1$x1) * rgamma(n, 2)
cases <- list(
  list(nm = "gamma errors, factor, qu 0.25, lsig -1", fo = y ~ x1 + x2 + f,
       qu = 0.25, lsig = -1, d = d1),
  list(nm = "gamma errors, factor, qu 0.95, lsig 0.5", fo = y ~ x1 + x2 + f,
       qu = 0.95, lsig = 0.5, d = d1))
for (cs in cases) {
  qf <- qgam(cs$fo, data = cs$d, qu = cs$qu, lsig = cs$lsig)
  m <- same_model(qf)
  r <- quiet(frm(bf(cs$fo, quantile = cs$qu, sigma = m$sigma),
                 family = quantile_elf(smoothing = m$smoothing), data = cs$d))
  ff <- r$v
  bq <- stats::coef(qf)
  bf_ <- fixef(ff)[, "Estimate"]
  cat(sprintf(paste0("  %s: smoothing %.5f | logLik qgam %.8f frm %.8f",
                     " (diff %.3g) | max |coef diff| / max(1, |coef|) %.3g |",
                     " conv %d warn %d %s\n"),
              cs$nm, m$smoothing, qgam_ll(qf, cs$d$y),
              as.numeric(logLik(ff)), qgam_ll(qf, cs$d$y) -
                as.numeric(logLik(ff)),
              max(abs(bq - bf_) / pmax(1, abs(bq))), ff$opt$convergence,
              length(r$w), paste(r$w, collapse = " | ")))
}

## 3b. one smooth case ----------------------------------------------------
cat("\n== 3b. qgam, smooth, sigma held\n")
set.seed(302)
d2 <- mgcv::gamSim(2, n = 500, dist = "normal", scale = 0.5,
                   verbose = FALSE)$data
d2 <- data.frame(x = d2$x, y = d2$y)
qf <- qgam(y ~ s(x, k = 15), data = d2, qu = 0.3, lsig = -1)
m <- same_model(qf)
r <- quiet(frm(bf(y ~ s(x, k = 15), quantile = 0.3, sigma = m$sigma),
               family = quantile_elf(smoothing = m$smoothing), data = d2,
               REML = TRUE))
ff <- r$v
muq <- as.numeric(stats::fitted(qf))
muf <- as.numeric(frm_linpred(ff, dpar = "mu"))
cat(sprintf(paste0("  gamSim(2) s(x, k = 15) qu 0.3 lsig -1: smoothing %.5f",
                   " | max |mu diff| / range %.3g | edf qgam %.3f | logLik",
                   " qgam %.6f frm %.6f | conv %d warn %d %s\n"),
            m$smoothing, max(abs(muq - muf)) / diff(range(muq)), sum(qf$edf),
            qgam_ll(qf, d2$y), as.numeric(logLik(ff)), ff$opt$convergence,
            length(r$w), paste(r$w, collapse = " | ")))

## 3c. rq() at small smoothing ---------------------------------------------
cat("\n== 3c. against rq(), in rq nid SE units\n")
for (s in 1:3) {
  set.seed(310 + s)
  n <- 500
  d <- data.frame(x = runif(n))
  d$y <- 1 + 2 * d$x + (0.3 + d$x) * rlnorm(n, 0, 0.75)
  for (tau in c(0.25, 0.9)) {
    rr <- rq(y ~ x, tau = tau, data = d)
    se <- summary(rr, se = "nid")$coefficients[, 2]
    out <- vapply(c(0.25, 0.05, 0.01, 0.002), function(lam) {
      r <- quiet(frm(bf(y ~ x, quantile = tau),
                     family = quantile_elf(smoothing = lam), data = d))
      if (length(r$w)) cat("   warn:", r$w, "\n")
      c(max(abs(fixef(r$v)[, "Estimate"] - coef(rr)) / se),
        r$v$opt$convergence)
    }, c(0, 0))
    cat(sprintf("  seed %d tau %.2f: max |coef - rq| / se at lambda 0.25 0.05 0.01 0.002: %s | conv %s\n",
                s, tau, paste(sprintf("%.3f", out[1, ]), collapse = " "),
                paste(out[2, ], collapse = " ")))
  }
}

## 3d. the bias formula on skewed residuals (population) ------------------
# maximize E log f(Y; mu, sigma) by integration, then compare mu* - q
# with -(pi^2/6) (lambda sigma*)^2 f'(q)/f(q) and with the exact smoothed
# quantile: P(Y - lambda sigma* L <= mu*) = p
cat("\n== 3d. population bias against the formula\n")
laws <- list(
  lognormal = list(d = function(y) dlnorm(y), q = function(p) qlnorm(p),
                   dl = function(y) -(1 + log(y)) / y, lo = 0, hi = Inf,
                   sd = sqrt((exp(1) - 1) * exp(1))),
  exponential = list(d = function(y) dexp(y), q = function(p) qexp(p),
                     dl = function(y) -1 + 0 * y, lo = 0, hi = Inf, sd = 1),
  gamma2 = list(d = function(y) dgamma(y, 2), q = function(p) qgamma(p, 2),
                dl = function(y) 1 / y - 1, lo = 0, hi = Inf,
                sd = sqrt(2)))
for (nm in names(laws)) {
  lw <- laws[[nm]]
  for (p in c(0.1, 0.5, 0.9)) {
    for (lam in c(0.25, 1)) {
      fam <- quantile_elf(smoothing = lam)
      q <- lw$q(p)
      el <- function(th) {
        mu <- th[1]
        sg <- exp(th[2])
        f <- function(y) {
          fam$lpdf(y, list(mu = mu, sigma = sg, quantile = p), list()) *
            lw$d(y)
        }
        br <- sort(unique(c(lw$lo, pmax(lw$lo, mu + c(-20, -1, 0, 1, 20) *
                                          lam * sg), lw$q(0.999))))
        v <- sum(vapply(seq_len(length(br) - 1), function(i) {
          stats::integrate(f, br[i], br[i + 1], rel.tol = 1e-12,
                           subdivisions = 1000)$value
        }, 0)) + stats::integrate(f, br[length(br)], Inf, rel.tol = 1e-12,
                                  subdivisions = 1000)$value
        -v
      }
      o <- stats::optim(c(q, log(lw$sd / 2)), el, method = "BFGS",
                        control = list(reltol = 1e-15, maxit = 500))
      o <- stats::optim(o$par, el, method = "Nelder-Mead",
                        control = list(reltol = 1e-15, maxit = 2000))
      mu <- o$par[1]
      s <- lam * exp(o$par[2])
      formula <- -(pi^2 / 6) * s^2 * lw$dl(q)
      # the exact smoothed quantile: solve E plogis((Y - m) / s) = 1 - p
      sm <- stats::uniroot(function(m) {
        stats::integrate(function(y) plogis((y - m) / s) * lw$d(y), lw$lo,
                         Inf, rel.tol = 1e-12, subdivisions = 1000)$value -
          (1 - p)
      }, q + c(-5, 5) * lw$sd, tol = 1e-12)$root
      cat(sprintf(paste0("  %-11s p %.1f lam %.2f: sigma* %.4f, mu* - q %+.5f",
                         " (%.4f sd) | formula %+.5f | smoothed quantile - q",
                         " %+.5f\n"),
                  nm, p, lam, exp(o$par[2]), mu - q, (mu - q) / lw$sd,
                  formula, sm - q))
    }
  }
}

## 3e. one large-sample fit, lognormal residuals ---------------------------
cat("\n== 3e. n = 40000 intercept-and-slope fit, lognormal residuals\n")
set.seed(320)
n <- 40000
d <- data.frame(x = runif(n))
d$y <- 1 + 2 * d$x + rlnorm(n)
for (p in c(0.5, 0.9)) {
  for (lam in c(0.25, 1)) {
    r <- quiet(frm(bf(y ~ x, quantile = p),
                   family = quantile_elf(smoothing = lam), data = d))
    b <- fixef(r$v)
    rr <- rq(y ~ x, tau = p, data = d)
    cat(sprintf(paste0("  p %.1f lam %.2f: intercept %.5f (se %.5f), truth",
                       " %.5f, rq %.5f | sigma %.4f | conv %d warn %d\n"),
                p, lam, b[1, "Estimate"], b[1, "Est.Error"], 1 + qlnorm(p),
                coef(rr)[1],
                as.numeric(asNamespace("frmtmb")$eval_dpars(r$v)[[1L]]$sigma[1]),
                r$v$opt$convergence, length(r$w)))
  }
}
cat("DONE\n")
