# Lane famend: one-step-ahead residuals on a censored count (3618) and
# with fractional weights on a continuous response (5726), before and
# after, against exact references.
#   Rscript dev/famend-osa.R <lib or "base"> [case ...]
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
cases <- a[-1]
want <- function(k) !length(cases) || k %in% cases
cap <- function(expr) {
  w <- character()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage")),
    error = function(e) structure(conditionMessage(e), class = "err"))
  list(v = v, w = w)
}
osa <- function(fit, ...) {
  r <- cap(residuals(fit, type = "osa", ...))
  for (x in r$w) cat("   W:", substr(gsub("\\s+", " ", x), 1, 150), "\n")
  if (inherits(r$v, "err")) {
    cat("   ERR:", substr(r$v, 1, 300), "\n")
    return(NULL)
  }
  r$v[, "Estimate"]
}
dp <- function(fit) lapply(frmtmb:::eval_dpars(fit)[[1L]], as.numeric)
# the exact randomized PIT of a count on the window (lo, hi], with the
# uniforms TMB draws (set.seed(123) over the uncensored rows)
win_ref <- function(y, Fq, pq, lo, hi, U) {
  Flo <- if (is.finite(lo)) Fq(lo) else 0
  Fhi <- if (is.finite(hi)) Fq(hi) else 1
  stats::qnorm((Fq(y) - Flo - U * pq(y)) / (Fhi - Flo))
}
rel <- function(r, ref) max(abs(r - ref)) / max(abs(ref))

if (want("pois")) {
  # 3618, the triage's construction (scratch-B/check1.R, 3544)
  set.seed(3)
  d <- data.frame(x = rnorm(100))
  d$y <- rpois(100, exp(1 + 0.3 * d$x))
  d$cn <- ifelse(d$y > 4, "right", "none")
  d$y[d$y > 4] <- 4
  fit <- frm(y | cens(cn) ~ x, family = poisson(), data = d)
  r <- osa(fit)
  if (!is.null(r)) {
    u <- d$cn == "none"
    m <- dp(fit)$mu
    U <- local({ set.seed(123); runif(sum(u)) })
    ref <- win_ref(d$y[u], function(q) ppois(q, m[u]),
                   function(q) dpois(q, m[u]), -Inf, 4, U)
    cat(sprintf(paste0("pois right at 4: censored rows NA %d of %d; ",
                       "uncensored finite %d of %d; max rel |r - exact| %.3g\n"),
                sum(is.na(r[!u])), sum(!u), sum(is.finite(r[u])), sum(u),
                rel(r[u], ref)))
  }
}
if (want("both")) {
  # left at 1 and right at 6, poisson, negbinomial and binomial
  set.seed(31)
  n <- 300
  d <- data.frame(x = rnorm(n))
  mu <- exp(1.2 + 0.4 * d$x)
  fams <- list(
    poisson = list(f = poisson(), y = rpois(n, mu),
                   F = function(q, p) ppois(q, p$mu),
                   P = function(q, p) dpois(q, p$mu)),
    negbinomial = list(f = negbinomial(), y = rnbinom(n, size = 3, mu = mu),
                       F = function(q, p) pnbinom(q, size = p$shape, mu = p$mu),
                       P = function(q, p) dnbinom(q, size = p$shape, mu = p$mu)),
    binomial = list(f = binomial(), y = rbinom(n, 10, plogis(-0.4 + d$x)),
                    F = function(q, p) pbinom(q, 10, p$mu),
                    P = function(q, p) dbinom(q, 10, p$mu)),
    zero_inflated_poisson = list(
      f = zero_inflated_poisson(),
      y = ifelse(runif(n) < 0.2, 0L, rpois(n, mu)),
      F = function(q, p) (q >= 0) * (p$zi + (1 - p$zi) * ppois(q, p$mu)),
      P = function(q, p) (q == 0) * p$zi + (1 - p$zi) * dpois(q, p$mu)))
  for (nm in names(fams)) {
    cs <- fams[[nm]]
    dd <- d
    dd$y <- pmin(pmax(cs$y, 1), 6)
    dd$cn <- ifelse(cs$y <= 1, "left", ifelse(cs$y > 6, "right", "none"))
    dd$k <- 10L
    fo <- if (nm == "binomial") bf(y | trials(k) + cens(cn) ~ x) else
      bf(y | cens(cn) ~ x)
    fit <- frm(fo, family = cs$f, data = dd)
    r <- osa(fit)
    if (is.null(r)) next
    u <- dd$cn == "none"
    p <- lapply(dp(fit), function(v) rep_len(v, n)[u])
    U <- local({ set.seed(123); runif(sum(u)) })
    ref <- win_ref(dd$y[u], function(q) cs$F(q, p), function(q) cs$P(q, p),
                   1, 6, U)
    inside <- r[u] >= qnorm(pmax(0, (cs$F(dd$y[u] - 1, p) - cs$F(1, p)) /
                                   (cs$F(6, p) - cs$F(1, p)))) - 1e-8 &
      r[u] <= qnorm((cs$F(dd$y[u], p) - cs$F(1, p)) /
                      (cs$F(6, p) - cs$F(1, p))) + 1e-8
    cat(sprintf(paste0("%s (1, 6]: left %d right %d; censored NA %d; ",
                       "max rel |r - exact| %.3g; in own step %d of %d\n"),
                nm, sum(dd$cn == "left"), sum(dd$cn == "right"),
                sum(is.na(r[!u])), rel(r[u], ref), sum(inside), sum(u)))
  }
}
if (want("re")) {
  # with (1 | g): the "cdf" route against oneStepGeneric's sum over the
  # support, and the row blocks against one call
  set.seed(41)
  G <- 30
  n <- 300
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = n / G)))
  yy <- rpois(n, exp(1 + 0.3 * d$x + rnorm(G, 0, 0.5)[d$g]))
  d$y <- pmin(pmax(yy, 1), 7)
  d$cn <- ifelse(yy <= 1, "left", ifelse(yy > 7, "right", "none"))
  fit <- frm(y | cens(cn) ~ x + (1 | g), family = poisson(), data = d)
  t0 <- proc.time()[[3]]
  r1 <- osa(fit)
  t1 <- proc.time()[[3]]
  r2 <- osa(fit, osa_method = "oneStepGeneric")
  t2 <- proc.time()[[3]]
  op <- options(frmtmb.osa_block = Inf)
  r3 <- osa(fit)
  t3 <- proc.time()[[3]]
  options(op)
  u <- d$cn == "none"
  if (is.null(r1) || is.null(r2) || is.null(r3)) {
    cat("re: refused\n")
  } else cat(sprintf(paste0("re poisson (1 | g), n %d, censored %d: cdf vs ",
                     "oneStepGeneric max |diff| %.3g median %.3g; blocks ",
                     "vs one call max |diff| %.3g; time cdf %.1fs generic ",
                     "%.1fs one call %.1fs\n"),
              n, sum(!u), max(abs(r1[u] - r2[u])),
              median(abs(r1[u] - r2[u])), max(abs(r1[u] - r3[u])),
              t1 - t0, t2 - t1, t3 - t2))
  if (!is.null(r1)) {
    cat("re: KS p of the cdf residuals against N(0, 1):",
        format(ks.test(r1[u], "pnorm")$p.value, digits = 3), "\n")
  }
}
if (want("cont")) {
  # a censored continuous response: "cdf" against the default
  # oneStepGeneric with exact quadrature, and against the exact PIT
  set.seed(51)
  n <- 200
  d <- data.frame(x = rnorm(n))
  y <- 1 + d$x + rnorm(n)
  d$y <- pmin(pmax(y, -0.5), 2.5)
  d$cn <- ifelse(y < -0.5, "left", ifelse(y > 2.5, "right", "none"))
  fit <- frm(y | cens(cn) ~ x, data = d)
  r1 <- osa(fit)
  r2 <- osa(fit, osa_method = "cdf")
  u <- d$cn == "none"
  p <- dp(fit)
  m <- rep_len(p$mu, n)[u]
  s <- rep_len(p$sigma, n)[u]
  ref <- qnorm((pnorm(d$y[u], m, s) - pnorm(-0.5, m, s)) /
                 (pnorm(2.5, m, s) - pnorm(-0.5, m, s)))
  cat(sprintf(paste0("cont gaussian [-0.5, 2.5]: generic vs exact %.3g, ",
                     "cdf vs exact %.3g (max rel)\n"),
              rel(r1[u], ref), rel(r2[u], ref)))
}
if (want("frac")) {
  # 5726: fractional weights on a continuous response (check1.R, 4939)
  set.seed(4)
  dw <- data.frame(x = rnorm(60), w = runif(60, 0.2, 1.8))
  dw$y <- 1 + dw$x + rnorm(60)
  fit <- frm(y | weights(w) ~ x, data = dw)
  r <- osa(fit)
  if (!is.null(r)) {
    p <- dp(fit)
    ref <- qnorm(pnorm(dw$y, rep_len(p$mu, 60), rep_len(p$sigma, 60)))
    f <- dw$w < 1
    cat(sprintf(paste0("frac gaussian: weight < 1 on %d rows, NA %d; ",
                       "max rel |r - qnorm(F(y))| fractional %.3g, ",
                       "others %.3g\n"), sum(f), sum(is.na(r)),
                if (all(is.finite(r[f]))) rel(r[f], ref[f]) else NA,
                rel(r[!f], ref[!f])))
  }
  # Gamma and a weight 0 row: NA there, and a count keeps NA with warning
  dw$yg <- rgamma(60, shape = 4, rate = 4 / exp(0.5 + 0.3 * dw$x))
  dw$w0 <- dw$w
  dw$w0[1:3] <- 0
  fg <- frm(yg | weights(w0) ~ x, family = Gamma("log"), data = dw)
  r <- osa(fg)
  if (!is.null(r)) {
    p <- dp(fg)
    sh <- rep_len(p$shape, 60)
    ref <- qnorm(pgamma(dw$yg, shape = sh, rate = sh / rep_len(p$mu, 60)))
    f <- dw$w0 > 0 & dw$w0 < 1
    cat(sprintf(paste0("frac Gamma: weight 0 rows NA %d of 3; fractional ",
                       "%d rows finite %d, max rel |r - exact| %.3g\n"),
                sum(is.na(r[1:3])), sum(f), sum(is.finite(r[f])),
                if (any(is.finite(r[f]))) rel(r[f], ref[f]) else NA))
  }
  dw$yc <- rpois(60, exp(0.5 + 0.3 * dw$x))
  fc <- frm(yc | weights(w) ~ x, family = poisson(), data = dw)
  r <- osa(fc)
  if (!is.null(r)) {
    cat(sprintf("frac poisson: NA on %d of %d fractional rows\n",
                sum(is.na(r[dw$w < 1])), sum(dw$w < 1)))
  }
  # under (1 | g): calibration of the fractional rows
  set.seed(61)
  G <- 40
  n <- 400
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), each = n / G)),
                  w = rep(c(1, 0.5), length.out = n))
  d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[d$g] + rnorm(n)
  fr <- frm(y | weights(w) ~ x + (1 | g), data = d)
  r <- osa(fr)
  if (!is.null(r)) {
    f <- d$w < 1
    cat(sprintf(paste0("frac gaussian (1 | g): fractional rows finite %d of ",
                       "%d, mean %.3f sd %.3f, KS p %.3g; weight-1 rows KS ",
                       "p %.3g\n"), sum(is.finite(r[f])), sum(f),
                mean(r[f]), sd(r[f]), ks.test(r[f], "pnorm")$p.value,
                ks.test(r[!f], "pnorm")$p.value))
  }
}
cat("\ndone\n")
