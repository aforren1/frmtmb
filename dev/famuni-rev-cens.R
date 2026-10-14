# Reviewer probe: cens() left/right/interval and trunc(lb, ub) on the
# three families, against log-likelihoods written by hand from
# independent CDFs, at the fit's own parameters.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
dp <- function(fit, nm) frm_linpred(fit, dpar = nm, type = "response")
report <- function(lab, fit, ll) {
  a <- as.numeric(logLik(fit))
  cat(sprintf("%-38s logLik %.10f hand %.10f rel diff %.2e conv %d\n", lab,
              a, ll, abs(a - ll) / abs(ll), fit$opt$convergence))
}
set.seed(1013)
n <- 400
x <- rnorm(n)

## frechet
nu <- 3
yt <- exp(1 + 0.3 * x) / gamma(1 - 1 / nu) * rexp(n)^(-1 / nu)
lb <- 0.8; ub <- 12
keep <- yt > lb & yt < ub
d <- data.frame(x = x[keep], yt = yt[keep])
m <- nrow(d)
d$cc <- sample(c(-1, 0, 1, 2), m, TRUE, prob = c(0.15, 0.55, 0.15, 0.15))
d$y <- d$yt
d$y2 <- d$yt
# interval rows: the bracket around the true value
d$y[d$cc == 2] <- pmax(lb, d$yt[d$cc == 2] * 0.8)
d$y2[d$cc == 2] <- pmin(ub, d$yt[d$cc == 2] * 1.3)
d$y[d$cc == -1] <- d$yt[d$cc == -1] * 1.2
d$y[d$cc == 1] <- d$yt[d$cc == 1] * 0.8
fit <- frm(y | cens(cc, y2) + trunc(lb = 0.8, ub = 12) ~ x,
           family = frechet(), data = d)
mu <- dp(fit, "mu"); nuh <- dp(fit, "nu")
s <- mu / gamma(1 - 1 / nuh)
Fr <- function(q) exp(-(q / s)^(-nuh))
lf <- log(nuh) - log(s) - (1 + nuh) * log(d$y / s) - (d$y / s)^(-nuh)
ll <- ifelse(d$cc == 0, lf,
      ifelse(d$cc == -1, log(Fr(d$y) - Fr(lb)),
      ifelse(d$cc == 1, log(Fr(ub) - Fr(d$y)),
             log(Fr(d$y2) - Fr(d$y))))) - log(Fr(ub) - Fr(lb))
report("frechet cens(-1,0,1,2) + trunc(lb, ub)", fit, sum(ll))
fit0 <- frm(y | cens(cc, y2) ~ x, family = frechet(), data = d)
mu <- dp(fit0, "mu"); nuh <- dp(fit0, "nu")
s <- mu / gamma(1 - 1 / nuh)
Fr <- function(q) exp(-(q / s)^(-nuh))
lf <- log(nuh) - log(s) - (1 + nuh) * log(d$y / s) - (d$y / s)^(-nuh)
ll <- ifelse(d$cc == 0, lf, ifelse(d$cc == -1, log(Fr(d$y)),
      ifelse(d$cc == 1, log1p(-Fr(d$y)), log(Fr(d$y2) - Fr(d$y)))))
report("frechet cens(-1,0,1,2)", fit0, sum(ll))

## GEV, with xi > 0 and with xi < 0
for (xi0 in c(0.2, -0.3)) {
  yt <- rg(n, 1 + 0.5 * x, 0.8, xi0)
  lb <- quantile(yt, 0.05); ub <- quantile(yt, 0.97)
  keep <- yt > lb & yt < ub
  d <- data.frame(x = x[keep], yt = yt[keep])
  m <- nrow(d)
  d$cc <- sample(c(-1, 0, 1, 2), m, TRUE, prob = c(0.15, 0.55, 0.15, 0.15))
  d$y <- d$yt; d$y2 <- d$yt
  d$y[d$cc == 2] <- pmax(lb, d$yt[d$cc == 2] - 0.3)
  d$y2[d$cc == 2] <- pmin(ub, d$yt[d$cc == 2] + 0.4)
  d$y[d$cc == -1] <- pmin(ub, d$yt[d$cc == -1] + 0.2)
  d$y[d$cc == 1] <- pmax(lb, d$yt[d$cc == 1] - 0.2)
  d$lb <- lb; d$ub <- ub
  fit <- frm(y | cens(cc, y2) + trunc(lb = lb, ub = ub) ~ x,
             family = gen_extreme_value(), data = d)
  mu <- dp(fit, "mu"); sg <- dp(fit, "sigma"); xi <- dp(fit, "xi")
  P <- function(q) evd::pgev(q, mu, sg, xi[1])
  lf <- evd::dgev(d$y, mu, sg, xi[1], log = TRUE)
  ll <- ifelse(d$cc == 0, lf,
        ifelse(d$cc == -1, log(P(d$y) - P(lb)),
        ifelse(d$cc == 1, log(P(ub) - P(d$y)),
               log(P(d$y2) - P(d$y))))) - log(P(ub) - P(lb))
  report(sprintf("GEV xi %.1f cens + trunc(lb, ub)", xi0), fit, sum(ll))
  # a trunc() bound outside the fitted support, where F is 0 or 1
  far <- if (xi0 > 0) min(d$y) - 50 else max(d$y2) + 50
  d$fb <- far
  f2 <- tryCatch(
    if (xi0 > 0) frm(y | cens(cc, y2) + trunc(lb = fb) ~ x,
                     family = gen_extreme_value(), data = d)
    else frm(y | cens(cc, y2) + trunc(ub = fb) ~ x,
             family = gen_extreme_value(), data = d),
    error = function(e) e)
  if (inherits(f2, "error")) {
    cat("GEV far bound: ERROR", conditionMessage(f2), "\n")
  } else {
    mu <- dp(f2, "mu"); sg <- dp(f2, "sigma"); xi <- dp(f2, "xi")
    P <- function(q) evd::pgev(q, mu, sg, xi[1])
    lf <- evd::dgev(d$y, mu, sg, xi[1], log = TRUE)
    ll <- ifelse(d$cc == 0, lf, ifelse(d$cc == -1, log(P(d$y)),
          ifelse(d$cc == 1, log1p(-P(d$y)), log(P(d$y2) - P(d$y)))))
    cat(sprintf("   far bound %.2f: end of support %.3f..%.3f\n", far,
                min(mu - sg / xi), max(mu - sg / xi)))
    report(sprintf("GEV xi %.1f cens + trunc at a far bound", xi0), f2,
           sum(ll))
  }
}

## discrete Weibull, frmtmb's inclusive convention
q <- plogis(1.2 + 0.4 * x)
k <- 1.1
yt <- pmax(ceiling((log1p(-runif(n)) / log(q))^(1 / k) - 1), 0)
lb <- 1; ub <- 15
keep <- yt >= lb & yt <= ub
d <- data.frame(x = x[keep], yt = yt[keep])
m <- nrow(d)
d$cc <- sample(c(-1, 0, 1, 2), m, TRUE, prob = c(0.15, 0.55, 0.15, 0.15))
d$y <- d$yt; d$y2 <- d$yt
d$y[d$cc == 2] <- pmax(lb, d$yt[d$cc == 2] - 1)
d$y2[d$cc == 2] <- pmin(ub, d$yt[d$cc == 2] + 2)
d$y[d$cc == -1] <- pmin(ub, d$yt[d$cc == -1] + 1)
d$y[d$cc == 1] <- pmax(lb, d$yt[d$cc == 1] - 1)
fit <- frm(y | cens(cc, y2) + trunc(lb = 1, ub = 15) ~ x,
           family = discrete_weibull(), data = d)
mu <- dp(fit, "mu"); kh <- dp(fit, "shape")
S <- function(v) { v <- rep_len(v, length(mu)); ifelse(v <= 0, 1, mu^(pmax(v, 0)^kh)) }  # P(Y >= v)
F <- function(v) 1 - S(v + 1)                            # P(Y <= v)
ll <- ifelse(d$cc == 0, log(S(d$y) - S(d$y + 1)),
      ifelse(d$cc == -1, log(F(d$y) - F(lb - 1)),
      ifelse(d$cc == 1, log(F(ub) - F(d$y - 1)),
             log(F(d$y2) - F(d$y - 1))))) - log(F(ub) - F(lb - 1))
report("dw cens(-1,0,1,2) + trunc(1, 15), incl.", fit, sum(ll))
llx <- ifelse(d$cc == 0, log(S(d$y) - S(d$y + 1)),
       ifelse(d$cc == -1, log(F(d$y) - F(lb - 1)),
       ifelse(d$cc == 1, log(F(ub) - F(d$y)),
              log(F(d$y2) - F(d$y))))) - log(F(ub) - F(lb - 1))
cat(sprintf("   brms's exclusive reading at the same parameters: %.6f (difference %.4f)\n",
            sum(llx), sum(llx) - sum(ll)))

cat("\n--- discrete Weibull, by censoring type\n")
for (code in c(0, -1, 1, 2)) {
  i <- d$cc == code
  sub <- d[i, ]
  cat(sprintf("code %2d rows %d: hand %.6f\n", code, sum(i), sum(ll[i])))
}
# one-type fits, each against its hand value
for (code in c(-1, 1, 2)) {
  dd <- d
  dd$cc[dd$cc != code & dd$cc != 0] <- 0
  dd$y[dd$cc == 0] <- dd$yt[dd$cc == 0]
  for (tr in c(FALSE, TRUE)) {
    f <- if (tr) frm(y | cens(cc, y2) + trunc(lb = 1, ub = 15) ~ x,
                     family = discrete_weibull(), data = dd)
         else frm(y | cens(cc, y2) ~ x, family = discrete_weibull(), data = dd)
    mu <- dp(f, "mu"); kh <- dp(f, "shape")
    S <- function(v) { v <- rep_len(v, length(mu)); ifelse(v <= 0, 1, mu^(pmax(v, 0)^kh)) }
    F <- function(v) 1 - S(v + 1)
    lo <- if (tr) F(lb - 1) else 0
    hi <- if (tr) F(ub) else 1
    l1 <- ifelse(dd$cc == 0, log(S(dd$y) - S(dd$y + 1)),
          ifelse(dd$cc == -1, log(F(dd$y) - lo),
          ifelse(dd$cc == 1, log(hi - F(dd$y - 1)),
                 log(F(dd$y2) - F(dd$y - 1))))) - log(hi - lo)
    # the alternative a left-censored row in a window might be scored by
    l2 <- ifelse(dd$cc == -1, log(F(dd$y) - (if (tr) F(lb) else 0)) -
                   log(hi - lo), l1)
    report(sprintf("dw code %d only, trunc %s (incl.)", code, tr), f, sum(l1))
    if (code == -1 && tr) {
      cat(sprintf("   with F(y) - F(lb) for the left rows: %.10f\n", sum(l2)))
    }
  }
}
## poisson, the same window, as a control on the core composition
q <- 3
yp <- rpois(n, exp(1 + 0.2 * x))
keep <- yp >= 1 & yp <= 15
dp2 <- data.frame(x = x[keep], y = yp[keep])
dp2$cc <- sample(c(-1, 0), nrow(dp2), TRUE, prob = c(0.3, 0.7))
fp <- frm(y | cens(cc) + trunc(lb = 1, ub = 15) ~ x, family = poisson(),
          data = dp2)
lam <- dp(fp, "mu")
Fp <- function(v) ppois(v, lam)
lpp <- ifelse(dp2$cc == 0, dpois(dp2$y, lam, log = TRUE),
              log(Fp(dp2$y) - Fp(0))) - log(Fp(15) - Fp(0))
lpp2 <- ifelse(dp2$cc == 0, dpois(dp2$y, lam, log = TRUE),
              log(Fp(dp2$y) - Fp(1))) - log(Fp(15) - Fp(0))
report("poisson left cens + trunc(1, 15), incl.", fp, sum(lpp))
cat(sprintf("   with F(y) - F(lb) for the left rows: %.10f\n", sum(lpp2)))
