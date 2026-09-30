# Reviewer, claims 3 to 5 on the lane build, own seeds and constructions:
#  A. the kappa warning: false-alarm side (interior data from an xbeta
#     with kappa > 0; zeros but no ones)
#  B. xbeta at high precision near its switch point: frmtmb's logLik
#     against a stats::pbeta() reference at the estimates, and the
#     reference optimum
#  C. zero_inflated_beta_binomial vs glmmTMB with (1 | g)
#  D. hurdle_cumulative: density vs a hand-written reference with disc
#     modeled (logit, the link brms cannot run), refusals and their
#     absent-condition twins, fitted()/conditional_effects() codes
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
capw <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(conditionMessage(e), class = "err")
    }),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  list(v = v, w = w)
}
ref_xbeta <- function(y, mu, phi, kappa) {
  a <- mu * phi; b <- (1 - mu) * phi; d <- 1 + 2 * kappa
  ifelse(y <= 0, stats::pbeta(kappa / d, a, b, log.p = TRUE),
         ifelse(y >= 1,
                stats::pbeta((1 + kappa) / d, a, b, lower.tail = FALSE,
                             log.p = TRUE),
                stats::dbeta((y + kappa) / d, a, b, log = TRUE) - log(d)))
}

cat("== A. kappa warning, false-alarm side ==\n")
# A1: xbeta truth with kappa = 1 and phi large enough that no row
# reaches either end: kappa is placed by the interior shape
set.seed(4101)
n <- 2000
z <- rbeta(n, 0.5 * 200, 0.5 * 200)
y <- 3 * z - 1
say("A1 data: min %.4f max %.4f, zeros %d, ones %d", min(y), max(y),
    sum(y <= 0), sum(y >= 1))
dA <- data.frame(y = y)
fA <- capw(frm(y ~ 1, family = xbeta(), data = dA))
fB <- frm(y ~ 1, family = Beta(), data = dA)
say("A1 warnings from frm(xbeta): %d%s", length(fA$w),
    if (length(fA$w)) paste0(" -> ", fA$w[1]) else "")
if (inherits(fA$v, "err")) { say("A1 frm(xbeta) ERROR: %s", substr(fA$v, 1, 120)); fA <- capw(frm(y ~ 1, family = xbeta(), data = dA, control = frmtmb_control(optimizer = "optim"))); say("A1 refit with optimizer = optim: warnings %d%s", length(fA$w), if (length(fA$w)) paste0(" -> ", substr(fA$w[1], 1, 90)) else "") }
ci <- confint(fA$v)
print(ci)
say("A1 kappa estimate %.4f (truth 1); logLik xbeta %.4f, Beta %.4f, diff %.2f",
    exp(ci["kappa_(Intercept)", "est"]), as.numeric(logLik(fA$v)),
    as.numeric(logLik(fB)), as.numeric(logLik(fA$v) - logLik(fB)))
# A2: many exact zeros, no ones
set.seed(4102)
z <- rbeta(n, 0.3 * 5, 0.7 * 5)
k <- 0.2
y <- pmin(pmax((1 + 2 * k) * z - k, 0), 1)
y[y >= 1] <- 0.999
say("A2 data: zeros %d, ones %d", sum(y == 0), sum(y == 1))
f2 <- capw(frm(y ~ 1, family = xbeta(), data = data.frame(y = y)))
say("A2 warnings: %d; kappa %.4f (truth 0.2)", length(f2$w),
    exp(confint(f2$v)["kappa_(Intercept)", "est"]))
# A3: kappa modeled by a predictor, with a response in (0, 1) only
set.seed(4103)
x <- rnorm(n)
z <- rbeta(n, 0.5 * 300, 0.5 * 300)
kk <- exp(0.2 * x)
y <- (1 + 2 * kk) * z - kk
say("A3 data (kappa ~ x): zeros %d ones %d", sum(y <= 0), sum(y >= 1))
f3 <- capw(frm(bf(y ~ 1, kappa ~ x), family = xbeta(), data = data.frame(y = y, x = x)))
say("A3 warnings: %d%s", length(f3$w),
    if (length(f3$w)) paste0(" -> ", substr(f3$w[1], 1, 90)) else "")
if (inherits(f3$v, "err")) say("A3 ERROR: %s", substr(f3$v, 1, 120)) else
  print(confint(f3$v))

cat("\n== B. xbeta at high precision, near the switch point ==\n")
for (phi in c(2e3, 2e4, 2e5)) {
  set.seed(4200)
  n <- 1000
  mu <- 0.047; kap <- 0.05
  z <- rbeta(n, mu * phi, (1 - mu) * phi)
  y <- pmin(pmax((1 + 2 * kap) * z - kap, 0), 1)
  dB <- data.frame(y = y)
  f <- capw(frm(y ~ 1, family = xbeta(), data = dB))
  if (inherits(f$v, "err")) { say("phi %g: ERROR %s", phi, f$v); next }
  est <- confint(f$v)[, "est"]
  m <- plogis(est[1]); ph <- exp(est[2]); ka <- exp(est[3])
  ref <- sum(ref_xbeta(y, m, ph, ka))
  ll <- as.numeric(logLik(f$v))
  # the reference's own optimum, started at frmtmb's
  nll <- function(p) -sum(ref_xbeta(y, plogis(p[1]), exp(p[2]), exp(p[3])))
  o <- optim(est, nll, method = "BFGS",
             control = list(reltol = 1e-14, maxit = 2000))
  o <- optim(o$par, nll, method = "Nelder-Mead",
             control = list(reltol = 1e-15, maxit = 5000))
  say(paste0("phi %g: zeros %d; est mu %.5f phi %.1f kappa %.5f; ",
             "logLik frmtmb %.8f ref-at-est %.8f diff %.2e; ",
             "ref optimum %.8f at mu %.5f phi %.1f kappa %.5f; warnings %d"),
      phi, sum(y == 0), m, ph, ka, ll, ref, ll - ref, -o$value,
      plogis(o$par[1]), exp(o$par[2]), exp(o$par[3]), length(f$w))
  if (length(f$w)) cat("   ", f$w, sep = "\n    ")
}

cat("\n== C. zero_inflated_beta_binomial vs glmmTMB, (1 | g) ==\n")
set.seed(4300)
n <- 1200
dC <- data.frame(x = rnorm(n), g = gl(60, n / 60),
                 tr = sample(3:20, n, TRUE))
u <- rnorm(60, 0, 0.5)[dC$g]
mu <- plogis(-0.4 + 0.6 * dC$x + u)
yb <- rbinom(n, dC$tr, rbeta(n, mu * 6, (1 - mu) * 6))
dC$y <- ifelse(runif(n) < plogis(-1.2 + 0.5 * dC$x), 0L, yb)
fC <- frm(bf(y | trials(tr) ~ x + (1 | g), zi ~ x),
          family = zero_inflated_beta_binomial(), data = dC)
gC <- glmmTMB::glmmTMB(cbind(y, tr - y) ~ x + (1 | g), ziformula = ~ x,
                       family = glmmTMB::betabinomial(), data = dC)
say("C logLik frmtmb %.10f glmmTMB %.10f diff %.2e",
    as.numeric(logLik(fC)), as.numeric(logLik(gC)),
    as.numeric(logLik(fC)) - as.numeric(logLik(gC)))
fe <- fixef(fC)
print(fe)
print(glmmTMB::fixef(gC))
say("C sd(g) frmtmb %.6f glmmTMB %.6f", sqrt(as.numeric(VarCorr(fC)[[1]]$sd^2 %||% NA)),
    attr(glmmTMB::VarCorr(gC)$cond$g, "stddev"))
say("C phi frmtmb %.6f glmmTMB %.6f",
    frmtmb:::eval_dpars(fC)[[1]]$phi[1], glmmTMB::sigma(gC))

cat("\n== D. hurdle_cumulative ==\n")
set.seed(4400)
n <- 1500
dD <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(30, n / 30))
disc <- exp(0.4 * dD$z)
u <- rlogis(n) / disc + 0.8 * dD$x
yc <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
dD$y <- ifelse(runif(n) < plogis(-0.6 + 0.5 * dD$x), 0L, yc)
say("D data: table %s", paste(table(dD$y), collapse = " "))
ref_hc <- function(y, eta, hu, disc, tau) {
  K <- length(tau) + 1
  Fm <- cbind(0, plogis(disc * (outer(rep(1, length(eta)), tau) - eta)), 1)
  p <- Fm[cbind(seq_along(y), pmax(y, 1) + 1)] - Fm[cbind(seq_along(y), pmax(y, 1))]
  ifelse(y == 0, log(hu), log1p(-hu) + log(p))
}
fD <- frm(bf(y ~ x, hu ~ x, disc ~ 0 + z), family = hurdle_cumulative(),
          data = dD)
dp <- frmtmb:::eval_dpars(fD)[[1]]
raw <- fD$estimates$tau_raw
tau <- cumsum(c(raw[1], exp(raw[-1])))
r <- sum(ref_hc(dD$y, dp$mu, dp$hu, dp$disc, tau))
say("D logit, disc ~ 0 + z: logLik %.10f reference %.10f rel diff %.2e",
    as.numeric(logLik(fD)), r, abs(as.numeric(logLik(fD)) - r) / abs(r))
say("D disc slope %.4f (truth 0.4); hu slope %.4f (truth 0.5)",
    fixef(fD)["disc_z", 1], fixef(fD)["hu_x", 1])
# fitted(): column 0 is hu, the rest (1 - hu) times cumulative probs
P <- fitted(fD)
say("D fitted dims %s; dimnames[[3]] %s", paste(dim(P), collapse = " x "),
    paste(dimnames(P)[[3]], collapse = ","))
P0 <- P[, "Estimate", ]
Fm <- cbind(0, plogis(dp$disc * (outer(rep(1, n), tau) - dp$mu)), 1)
Pref <- cbind(dp$hu, (1 - dp$hu) * (Fm[, -1] - Fm[, -ncol(Fm)]))
say("D fitted vs reference: max abs diff %.2e", max(abs(P0 - Pref)))
ce <- conditional_effects(fD, effects = "x", categorical = FALSE)[[1]]
nd <- ce[, c("x", "z")]
Pn <- fitted(fD, newdata = nd, re_formula = NA)[, "Estimate", ]
say("D ce(categorical = FALSE) vs sum_k k P_k with codes 0..4: max diff %.2e; vs codes 1..5 %.2e",
    max(abs(ce$estimate__ - as.vector(Pn %*% (0:4)))),
    max(abs(ce$estimate__ - as.vector(Pn %*% (1:5)))))
# refusals and their absent twins
dD$grp <- factor(ifelse(dD$z > 0, "p", "n"))
dD$fc <- factor(sample(c("a", "b"), n, TRUE))
chk <- function(lab, expr) {
  r <- capw(expr)
  say("  %-42s %s%s", lab,
      if (inherits(r$v, "err")) paste("REFUSED:", substr(r$v, 1, 110)) else
        paste("fits, logLik", format(as.numeric(logLik(r$v)), digits = 10)),
      if (length(r$w)) paste(" [warn:", substr(r$w[1], 1, 60), "]") else "")
}
chk("thres(gr = grp)", frm(y | thres(gr = grp) ~ x, data = dD,
                           family = hurdle_cumulative()))
chk("thres(4) (absent twin)", frm(y | thres(4) ~ x, data = dD,
                                  family = hurdle_cumulative(),
                                  prior = set_prior("normal(0, 3)",
                                                    class = "Intercept")))
chk("cs(fc)", frm(y ~ x + cs(fc), data = dD, family = hurdle_cumulative()))
chk("cse(fc)", frm(y ~ x + cse(fc), data = dD, family = hurdle_cumulative()))
chk("fc (absent twin)", frm(y ~ x + fc, data = dD,
                            family = hurdle_cumulative()))
chk("threshold = 'equidistant'", frm(y ~ x, data = dD,
                    family = hurdle_cumulative(threshold = "equidistant")))
chk("threshold = 'flexible' (absent twin)", frm(y ~ x, data = dD,
                    family = hurdle_cumulative(threshold = "flexible")))
chk("disc ~ 1 + z (intercept, unidentified)",
    frm(bf(y ~ x, disc ~ 1 + z), data = dD, family = hurdle_cumulative()))
chk("disc ~ 0 + z, (1 | g) on mu",
    frm(bf(y ~ x + (1 | g), disc ~ 0 + z), data = dD,
        family = hurdle_cumulative()))
chk("hu ~ x + (1 | g)",
    frm(bf(y ~ x, hu ~ x + (1 | g)), data = dD, family = hurdle_cumulative()))
chk("mixture(hurdle_cumulative, hurdle_cumulative)",
    frm(y ~ x, data = dD, family = mixture(hurdle_cumulative(),
                                           hurdle_cumulative())))
chk("response 0..1 only", frm(y ~ x, data = transform(dD, y = pmin(y, 1)),
                               family = hurdle_cumulative()))
chk("response 1..4 (no zeros)", frm(y ~ x, data = transform(dD, y = pmax(y, 1)),
                                    family = hurdle_cumulative()))
dD$yo <- factor(c("none", "low", "mid", "high", "top")[dD$y + 1],
                levels = c("none", "low", "mid", "high", "top"), ordered = TRUE)
fo <- frm(bf(yo ~ x, hu ~ x, disc ~ 0 + z), family = hurdle_cumulative(),
          data = dD)
say("D ordered factor response: logLik %.10f vs integer %.10f",
    as.numeric(logLik(fo)), as.numeric(logLik(fD)))
s <- simulate(fo, nsim = 1, seed = 1)[[1]]
say("D simulate() on the factor fit: class %s, levels %s, table %s",
    paste(class(s), collapse = "/"), paste(levels(s), collapse = ","),
    paste(table(s), collapse = " "))
say("D predict() column names: %s",
    paste(colnames(predict(fo))[1:5], collapse = ", "))
