# Review lane osa: residuals(type = "osa") against the EXACT one-step
# predictive, on designs and seeds of the reviewer's own.
#
#   Rscript dev/review-osa-exact.R <seed> [lib]
#
# Design: 20 groups of 10 rows in time order, x ~ N(0, 1),
# eta = 0.2 + 0.5 x + e. Under "re", e is a group intercept, sd 1.0
# (fit (1 | g)). Under "ar", e is a latent AR(1) per group, coefficient
# 0.5, innovation sd 0.8 (fit ar(t, g, cov = TRUE)). Trials vary by row
# (1 to 12) for binomial and beta_binomial.
#
# Exact: at the fit's own estimates, the one-step predictive of row i is
# the integral of F(y_i | e) against p(e | earlier rows of its group).
# "re": a grid of 1601 points over +-8 sd. "ar": a forward filter on a grid
# of 801 states over +-8 stationary sd. Every density and CDF here is
# written from stats:: functions, with no frmtmb code; the zero-inflated
# and hurdle forms are brms's. With TMB's own uniforms
# (set.seed(123); runif(n)) the exact residual is qnorm(F(y) - U p(y)).
args <- commandArgs(TRUE)
seed <- as.integer(args[1])
lib1 <- if (length(args) >= 2) args[2] else "C:/Users/adf44/source/r/wt-osa-lib"
.libPaths(c(if (lib1 != "none") lib1, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")

G <- 20L
K <- 10L
n <- G * K
U <- local({ set.seed(123); runif(n) })

lbb <- function(y, k, a, b) {
  lchoose(k, y) + lbeta(y + a, k - y + b) - lbeta(a, b)
}
# p(y | eta) and F(y | eta) for one row over a vector of eta
fams <- list(
  bernoulli = list(
    fam = bernoulli(), inv = plogis,
    p = function(y, m, dp, k) dbinom(y, 1, m),
    F = function(y, m, dp, k) pbinom(y, 1, m),
    r = function(m, dp, k) rbinom(length(m), 1, m)),
  binomial = list(
    fam = binomial(), inv = plogis,
    p = function(y, m, dp, k) dbinom(y, k, m),
    F = function(y, m, dp, k) pbinom(y, k, m),
    r = function(m, dp, k) rbinom(length(m), k, m)),
  poisson = list(
    fam = poisson(), inv = exp,
    p = function(y, m, dp, k) dpois(y, m),
    F = function(y, m, dp, k) ppois(y, m),
    r = function(m, dp, k) rpois(length(m), m)),
  negbinomial = list(
    fam = negbinomial(), inv = exp,
    p = function(y, m, dp, k) dnbinom(y, size = dp$shape, mu = m),
    F = function(y, m, dp, k) pnbinom(y, size = dp$shape, mu = m),
    r = function(m, dp, k) rnbinom(length(m), size = 3, mu = m)),
  beta_binomial = list(
    fam = beta_binomial(), inv = plogis,
    p = function(y, m, dp, k) exp(lbb(y, k, m * dp$phi, (1 - m) * dp$phi)),
    F = function(y, m, dp, k) {
      if (y < 0) return(0 * m)
      Reduce(`+`, lapply(0:y, function(j) {
        exp(lbb(j, k, m * dp$phi, (1 - m) * dp$phi))
      }))
    },
    r = function(m, dp, k) rbinom(length(m), k, rbeta(length(m), m * 8,
                                                       (1 - m) * 8))),
  zero_inflated_poisson = list(
    fam = zero_inflated_poisson(), inv = exp,
    p = function(y, m, dp, k) dp$zi * (y == 0) + (1 - dp$zi) * dpois(y, m),
    F = function(y, m, dp, k) {
      if (y < 0) return(0 * m)
      dp$zi + (1 - dp$zi) * ppois(y, m)
    },
    r = function(m, dp, k) ifelse(runif(length(m)) < 0.3, 0,
                                  rpois(length(m), m))),
  hurdle_poisson = list(
    fam = hurdle_poisson(), inv = exp,
    p = function(y, m, dp, k) {
      if (y == 0) return(dp$hu + 0 * m)
      (1 - dp$hu) * dpois(y, m) / (1 - dpois(0, m))
    },
    F = function(y, m, dp, k) {
      if (y < 0) return(0 * m)
      dp$hu + (1 - dp$hu) * (ppois(y, m) - dpois(0, m)) / (1 - dpois(0, m))
    },
    r = function(m, dp, k) {
      p0 <- dpois(0, m)
      pos <- qpois(runif(length(m), p0, 1), m)
      ifelse(runif(length(m)) < 0.3, 0, pmax(pos, 1))
    }))

make_data <- function(fam, design) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(G), each = K)), t = rep(1:K, G),
                  x = rnorm(n), k = sample(1:12, n, replace = TRUE))
  e <- if (design == "re") {
    rnorm(G, 0, 1.0)[d$g]
  } else {
    as.vector(replicate(G, arima.sim(list(ar = 0.5), K, sd = 0.8)))
  }
  m <- fams[[fam]]$inv(0.2 + 0.5 * d$x + e)
  d$y <- fams[[fam]]$r(m, NULL, d$k)
  d
}

form_of <- function(fam, design) {
  rhs <- if (design == "re") "x + (1 | g)" else "x + ar(t, g, cov = TRUE)"
  lhs <- if (fam %in% c("binomial", "beta_binomial")) "y | trials(k)" else "y"
  stats::as.formula(paste(lhs, "~", rhs))
}

exact_re <- function(fit, d, ff, dp) {
  eta <- fit$estimates$beta[["(Intercept)"]] + fit$estimates$beta[["x"]] * d$x
  sdu <- exp(fit$estimates$theta[[1]])
  u <- seq(-8, 8, length.out = 1601) * sdu
  w0 <- dnorm(u, 0, sdu)
  lo <- hi <- numeric(n)
  for (gi in split(seq_len(n), d$g)) {
    w <- w0 / sum(w0)
    for (i in gi) {
      m <- ff$inv(eta[i] + u)
      p <- ff$p(d$y[i], m, dp, d$k[i])
      Fh <- ff$F(d$y[i], m, dp, d$k[i])
      hi[i] <- sum(w * Fh)
      lo[i] <- sum(w * (Fh - p))
      w <- w * p
      w <- w / sum(w)
    }
  }
  list(lo = lo, hi = hi)
}

exact_ar <- function(fit, d, ff, dp) {
  eta <- fit$estimates$beta[["(Intercept)"]] + fit$estimates$beta[["x"]] * d$x
  th <- fit$estimates$theta
  sderr <- exp(th[[1]])
  phi <- frmtmb:::autocor_pacf(th[[2]])
  sds <- 1 / sqrt(1 - phi^2)
  s <- seq(-8, 8, length.out = 801) * sds
  Kt <- outer(s, s, function(a, b) dnorm(b, phi * a, 1))
  w0 <- dnorm(s, 0, sds)
  lo <- hi <- numeric(n)
  for (gi in split(seq_len(n), d$g)) {
    w <- w0 / sum(w0)
    for (i in gi) {
      m <- ff$inv(eta[i] + sderr * s)
      p <- ff$p(d$y[i], m, dp, d$k[i])
      Fh <- ff$F(d$y[i], m, dp, d$k[i])
      hi[i] <- sum(w * Fh)
      lo[i] <- sum(w * (Fh - p))
      w <- as.vector((w * p) %*% Kt)
      w <- w / sum(w)
    }
  }
  list(lo = lo, hi = hi, phi = phi, sderr = sderr)
}

cat("BEGIN review-osa-exact seed", seed, "\n")
for (design in c("re", "ar")) {
  for (fam in names(fams)) {
    ff <- fams[[fam]]
    d <- make_data(fam, design)
    fit <- tryCatch(suppressWarnings(frm(form_of(fam, design), family = ff$fam,
                                         data = d)), error = function(e) e)
    if (inherits(fit, "error")) {
      cat(sprintf("%-3s %-22s FIT ERROR %s\n", design, fam,
                  conditionMessage(fit)))
      next
    }
    dpe <- lapply(frmtmb:::eval_dpars(fit)[[1]], as.numeric)
    dp <- list(shape = dpe$shape[1], phi = dpe$phi[1], zi = dpe$zi[1],
               hu = dpe$hu[1])
    t0 <- proc.time()[["elapsed"]]
    wn <- character(0)
    r <- tryCatch(withCallingHandlers(
      residuals(fit, type = "osa")[, "Estimate"],
      warning = function(w) {
        wn <<- c(wn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    if (inherits(r, "error")) {
      cat(sprintf("%-3s %-22s OSA ERROR %s\n", design, fam,
                  substr(conditionMessage(r), 1, 200)))
      next
    }
    ex <- if (design == "re") exact_re(fit, d, ff, dp) else
      exact_ar(fit, d, ff, dp)
    ez <- qnorm(ex$hi - U * (ex$hi - ex$lo))
    df <- abs(r - ez)
    ok <- is.finite(df)
    j <- which.max(replace(df, !ok, -1))
    extra <- if (design == "ar") {
      sprintf(" ar %.2f sderr %.2f", ex$phi, ex$sderr)
    } else sprintf(" sd %.2f", exp(fit$estimates$theta[[1]]))
    cat(sprintf(paste0("%-3s %-22s nan %d | max %.3g med %.3g p95 %.3g ",
                       ">0.01 %d >0.05 %d | at row %d t %d y %g k %d ",
                       "r %.3f exact %.3f step [%.3f, %.3f] | mean %.3f ",
                       "sd %.3f ks.p %.3g | exact ks.p %.3g |%s | %.1fs%s\n"),
                design, fam, sum(!is.finite(r)), max(df[ok]),
                median(df[ok]), quantile(df[ok], 0.95), sum(df[ok] > 0.01),
                sum(df[ok] > 0.05), j, d$t[j], d$y[j], d$k[j], r[j], ez[j],
                qnorm(ex$lo[j]), qnorm(ex$hi[j]), mean(r[ok]), sd(r[ok]),
                suppressWarnings(ks.test(r[ok], "pnorm")$p.value),
                suppressWarnings(ks.test(ez[is.finite(ez)], "pnorm")$p.value),
                extra, el,
                if (length(wn)) paste0(" WARN: ", substr(wn[1], 1, 120)) else ""))
  }
}
cat("END review-osa-exact\n")
