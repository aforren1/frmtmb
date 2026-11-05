# Fix round 2, item 3: the density's join, the CDF, the discrete check. Build lib2.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
cat("--- the join: lgamma_count_shift against the exact sum, value and gradient, s = 23 to 26\n")
ex <- function(a, s) sum(log(a + seq_len(s) - 1)); exg <- function(a, s) sum(1 / (a + seq_len(s) - 1))
worst <- 0; worstg <- 0; rows <- NULL
for (s in c(1, 23, 24, 25, 26, 100)) for (a in c(1e-8, 0.01, 1, 100, 1e6, 1e12, 1e30)) {
  v <- as.numeric(ns$lgamma_count_shift(a, s))
  tp <- RTMB::MakeTape(function(p) ns$lgamma_count_shift(p, s), a)
  vt <- as.numeric(tp(a)); g <- as.numeric(tp$jacobian(a))
  rv <- abs(v - ex(a, s)) / max(1, abs(ex(a, s))); rt <- abs(vt - ex(a, s)) / max(1, abs(ex(a, s))); rg <- abs(g - exg(a, s)) / abs(exg(a, s))
  rows <- rbind(rows, data.frame(s = s, a = a, val = max(rv, rt), grad = rg))
}
print(aggregate(cbind(val, grad) ~ s, rows, max))
print(rows[which.max(rows$grad), ]); print(rows[rows$s %in% 24:25 & rows$a %in% c(1e-8, 1e30), ])
cat("--- rows with 24 and 25 trials in one fit, against dbetabinom at moderate phi and the finite product at large phi\n")
fp <- function(y, n, a, b) { phi <- a + b; lchoose(n, y) + (if (y > 0) sum(log(a + 0:(y - 1))) else 0) + (if (n - y > 0) sum(log(b + 0:(n - y - 1))) else 0) - sum(log(phi + 0:(n - 1))) }
set.seed(13001); nn <- 200; yy <- rbinom(nn, rep(c(24, 25), length.out = nn), 0.4); sz <- rep(c(24, 25), length.out = nn)
for (lp in c(-3, 2, 15, 39, 80)) {
  phi <- exp(lp); a <- 0.4 * phi; b <- 0.6 * phi
  mine <- ns$beta_binomial_lpmf(yy, sz, rep(a, nn), rep(b, nn), rep(phi, nn))
  h <- mapply(fp, yy, sz, MoreArgs = list(a = a, b = b))
  cat(sprintf("log phi %3d: max |lpmf - product| %.2e, sum %.6f\n", lp, max(abs(mine - h)), sum(mine)))
}
cat("--- the CDF against the exact sum, doubles and tape; edges; per-row trials\n")
excdf <- function(q, n, a, b) { if (q < 0) return(-Inf); if (q >= n) return(0); log(sum(exp(sapply(0:q, fp, n = n, a = a, b = b)))) }
bad <- 0L; maxd <- 0; maxg <- 0; cases <- 0L
for (n in c(1, 5, 10, 24, 25, 60)) for (lp in c(-3, 0, 2, 15, 30, 39, 60, 100)) for (m in c(0.02, 0.5, 0.98)) {
  phi <- exp(lp); qs <- sort(unique(c(-1, 0, 1, floor(n / 2), n - 1, n, n + 1))); N <- length(qs)
  cd <- ns$beta_binomial_cdf_log(qs, rep(n, N), rep(m * phi, N), rep((1 - m) * phi, N)); v <- as.numeric(cd$lF)
  hs <- mapply(function(q) { if (q < 0) return(0); if (q >= n) return(-Inf); log(sum(exp(sapply((q + 1):n, fp, n = n, a = m * phi, b = (1 - m) * phi)))) }, qs); oks <- is.finite(hs) & is.finite(as.numeric(cd$lS)); maxs <- max(c(if (exists("maxs")) maxs, abs(as.numeric(cd$lS)[oks] - hs[oks]) / pmax(1, abs(hs[oks]))), 0); if (any(is.finite(hs) != is.finite(as.numeric(cd$lS)))) bad <- bad + 1L
  h <- mapply(excdf, qs, MoreArgs = list(n = n, a = m * phi, b = (1 - m) * phi))
  ok <- is.finite(h) & is.finite(v)
  d <- max(abs(v[ok] - h[ok]) / pmax(1e-300, pmax(1, abs(h[ok]))), 0)
  if (any(is.finite(h) != is.finite(v))) { bad <- bad + 1L; if (bad <= 4) cat(" finite mismatch n", n, "lp", lp, "m", m, ":", v, "|", h, "\n") }
  maxd <- max(maxd, d); cases <- cases + 1L
  # the derivative in logit mu on the tape, against differences of the exact sum
  q0 <- floor(n / 2)
  if (q0 < n && lp <= 60) {
    tp <- RTMB::MakeTape(function(p) { pr <- RTMB::plogis(p); ns$beta_binomial_cdf_log(q0, n, pr * phi, (1 - pr) * phi)$lF }, qlogis(m))
    g <- as.numeric(tp$jacobian(qlogis(m))); e <- 1e-5
    fd <- (excdf(q0, n, plogis(qlogis(m) + e) * phi, (1 - plogis(qlogis(m) + e)) * phi) - excdf(q0, n, plogis(qlogis(m) - e) * phi, (1 - plogis(qlogis(m) - e)) * phi)) / (2 * e)
    if (is.finite(fd) && abs(fd) > 1e-6) maxg <- max(maxg, abs(g - fd) / abs(fd))
  }
}
cat(sprintf("cases %d, largest relative difference %.2e, finite mismatches %d, largest relative gradient difference %.2e\n", cases, maxd, bad, maxg))
qm <- c(3, 12, 30); szm <- c(10, 24, 40); vm <- as.numeric(ns$beta_binomial_cdf_log(qm, szm, 0.3 * c(5, 5, 5), 0.7 * c(5, 5, 5))$lF)
cat("per-row trials:", max(abs(vm - mapply(excdf, qm, szm, MoreArgs = list(a = 1.5, b = 3.5)))), "\n")
cat("--- censored and truncated fits, no overdispersion, against the hand likelihood at the estimate\n")
for (s in c(9, 6)) {
  set.seed(s); n <- 300
  d <- data.frame(x = round(rnorm(n), 3), ntr = 10L); d$y <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x))
  d$ce <- ifelse(d$y >= 8, "right", ifelse(d$y <= 1, "left", "none")); d$yc <- pmin(pmax(d$y, 1), 8)
  dtr <- d[d$y <= 8, ]
  r1 <- tt(frm(bf(yc | trials(ntr) + cens(ce) ~ x), family = beta_binomial(), data = d))
  r2 <- tt(frm(bf(y | trials(ntr) + trunc(ub = 8) ~ x), family = beta_binomial(), data = dtr))
  r3 <- tt(frm(bf(y | trials(ntr) + trunc(lb = 1, ub = 8) ~ x), family = beta_binomial(), data = dtr[dtr$y >= 1, ]))
  r4 <- tt(frm(bf(yc | trials(ntr) + cens(ce) ~ x), family = zero_inflated_beta_binomial(), data = d))
  hand <- function(r, type) {
    b <- r$estimates$beta; phi <- exp(r$estimates$betad[1])
    dd <- switch(type, cens = d, ub = dtr, lub = dtr[dtr$y >= 1, ])
    p <- plogis(b[1] + b[2] * dd$x)
    sum(vapply(seq_len(nrow(dd)), function(i) {
      a <- p[i] * phi; bb <- (1 - p[i]) * phi
      switch(type,
        cens = if (dd$ce[i] == "none") fp(dd$yc[i], 10, a, bb) else if (dd$ce[i] == "right") log(sum(exp(sapply(dd$yc[i]:10, fp, n = 10, a = a, b = bb)))) else log(sum(exp(sapply(0:dd$yc[i], fp, n = 10, a = a, b = bb)))),
        ub = fp(dd$y[i], 10, a, bb) - excdf(8, 10, a, bb),
        lub = fp(dd$y[i], 10, a, bb) - log(sum(exp(sapply(1:8, fp, n = 10, a = a, b = bb)))))
    }, 0))
  }
  for (x in list(list("cens right and left", r1, "cens"), list("trunc ub", r2, "ub"), list("trunc lb ub", r3, "lub"))) {
    r <- x[[2]]
    if (iserr(r)) cat(sprintf("seed %d %-20s ERROR %s\n", s, x[[1]], substr(r$msg, 1, 150))) else
      cat(sprintf("seed %d %-20s ll %.6f code %d log phi %.2f, hand %.6f, diff %.2e\n", s, x[[1]], as.numeric(logLik(r)), r$opt$convergence, r$estimates$betad[1], hand(r, x[[3]]), as.numeric(logLik(r)) - hand(r, x[[3]])))
  }
  cat(sprintf("seed %d zi cens: %s\n", s, if (iserr(r4)) substr(r4$msg, 1, 150) else sprintf("ll %.6f code %d log phi %.2f", as.numeric(logLik(r4)), r4$opt$convergence, r4$estimates$betad[1])))
}
cat("--- memory and time of the table at many trials: n = 2000, trials = 2000, censored rows only\n")
set.seed(13002); n <- 2000; dl <- data.frame(x = rnorm(n), ntr = 2000L); dl$y <- rbinom(n, 2000, plogis(0.3 + 0.2 * dl$x)); dl$ce <- ifelse(runif(n) < 0.05, "right", "none")
t1 <- system.time(rl <- tt(frm(bf(y | trials(ntr) + cens(ce) ~ x), family = beta_binomial(), data = dl)))[["elapsed"]]
cat(sprintf("%s: %.1f s, %s, peak via gc %.0f MB\n", arm, t1, if (iserr(rl)) substr(rl$msg, 1, 120) else format(as.numeric(logLik(rl)), digits = 10), sum(gc()[, 6])))
if (arm == "before") quit("no")
cat("\n--- the discrete check: try to make it fire falsely\n")
fires <- function(r) any(grepl("every response is discrete", attr(r, "w"), fixed = TRUE))
sh <- function(lab, r) cat(sprintf("%-38s %s fires=%s\n", lab, if (iserr(r)) paste("ERR", substr(r$msg, 1, 110)) else format(as.numeric(logLik(r)), digits = 8), if (iserr(r)) NA else fires(r)))
suppressPackageStartupMessages({ library(frmtmb.learn); library(frmtmb.eam); library(frmtmb.latent) })
set.seed(13003)
dl <- data.frame(id = rep(1:6, each = 40), trial = rep(1:40, times = 6), rt = 0.4 + rnorm(240, 0, 0.004), choice = rbinom(240, 1, 0.5), pay1 = rbinom(240, 1, 0.7), pay2 = rbinom(240, 1, 0.3))
sh("rlddm, tight rt (rt is a density)", tt(frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1), family = rlddm(subject = id, trial = trial), data = dl)))
dl$rt <- 0.3 + rlnorm(240, -4, 0.2)
sh("rlddm, short rt", tt(frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1), family = rlddm(subject = id, trial = trial), data = dl)))
print(frm_learn_families())
set.seed(13004); n <- 120
dz <- data.frame(g = factor(rep(1:12, each = 10)), y0 = 0L, x = rnorm(n))
sh("ZIP, y all 0, RE on mu (Laplace)", tt(frm(y0 ~ 1 + (1 | g), family = zero_inflated_poisson(), data = dz)))
sh("hurdle_poisson y all 0, RE", tt(frm(bf(y0 ~ 1 + (1 | g), hu ~ 1), family = hurdle_poisson(), data = dz)))
dz$y1 <- 1L
sh("bernoulli y all 1, RE", tt(frm(y1 ~ 1 + (1 | g), family = bernoulli(), data = dz)))
sh("mixture of poissons, y all 0, RE", tt(frm(y0 ~ 1 + (1 | g), family = mixture(poisson(), poisson()), data = dz)))
dz$ys <- rbinom(n, 1, 0.5 * (dz$x > 0) + 0.02); dz$y2 <- as.integer(dz$x > 0)
sh("bernoulli separation", tt(frm(y2 ~ x, family = bernoulli(), data = dz)))
sh("bernoulli separation, RE", tt(frm(y2 ~ x + (1 | g), family = bernoulli(), data = dz)))
sh("poisson weights 5", tt(frm(y0 | weights(rep(5, n)) ~ 1, family = poisson(), data = dz)))
sh("cumulative separation (not 'discrete')", tt(frm(ordered(y2) ~ x, family = cumulative(), data = dz)))
sh("categorical separation", tt(frm(factor(y2 + 1L + as.integer(x > 1)) ~ x, family = categorical(), data = dz)))
sh("trunc: poisson, ub 3, y = 0..3", tt(frm(yp | trunc(ub = 3) ~ 1, family = poisson(), data = transform(dz, yp = pmin(rpois(n, 1), 3)))))
sh("binomial trials 1, p exact", tt(frm(y2 | trials(rep(1L, n)) ~ x, family = binomial(), data = dz)))
sh("two discrete responses", tt(frm(mvbf(bf(y2 ~ x, family = bernoulli()), bf(y0 ~ 1, family = poisson()), rescor = FALSE), data = dz)))
sh("quadrature, ZIP y all 0", tt(frm(y0 ~ 1 + (1 | g), family = zero_inflated_poisson(), data = dz, quadrature = TRUE)))
ydc <- rep(0:1, each = 60); sh("discrete_weibull perfect", tt(frm(yy ~ 1, family = discrete_weibull(), data = data.frame(yy = rep(0L, 60)))))
sh("lca", tt({ set.seed(1); yl <- rbinom(60, 1, 0.5); frm(y ~ 1, family = frmtmb.latent::lca(2, bernoulli()), data = data.frame(y = yl, id = factor(rep(1:20, each = 3)))) }))
cat("check fires on what it is for: a raw density as a custom family, seed 9\n")
raw <- frmtmb_family("rawbb", dpars = c("mu", "phi"), links = list(mu = "logit", phi = "log"), type = "discrete", required_aterms = "trials",
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) RTMBdist::dbetabinom(y, ns$aterm_trials(aterms, "rawbb"), dpars$mu * dpars$phi, (1 - dpars$mu) * dpars$phi, log = TRUE),
  init_dpars = list(mu = function(y, aterms) 0.5, phi = function(y, aterms) 5))
set.seed(9); d9 <- data.frame(x = round(rnorm(300), 3), ntr = 10L); d9$y <- rbinom(300, 10, plogis(0.3 + 0.5 * d9$x))
sh("raw density, seed 9", tt(frm(bf(y | trials(ntr) ~ x), family = raw, data = d9)))
