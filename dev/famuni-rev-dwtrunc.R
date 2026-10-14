# Reviewer probe: discrete_weibull() with trunc(), no cens(), against a
# hand log-likelihood at the fit's own parameters.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
dp <- function(fit, nm) {
  e <- frm_linpred(fit, dpar = nm, type = "link")
  if (nm == "mu") plogis(e) else exp(e)
}
dp_resp <- function(fit, nm) frm_linpred(fit, dpar = nm, type = "response")
set.seed(1014)
n <- 400
x <- rnorm(n)
q <- plogis(1.2 + 0.4 * x)
yt <- pmax(ceiling((log1p(-runif(n)) / log(q))^(1 / 1.1) - 1), 0)
for (cfg in list(c(1, Inf), c(0, 15), c(1, 15), c(2, 9))) {
  lb <- cfg[1]; ub <- cfg[2]
  keep <- yt >= lb & yt <= ub
  d <- data.frame(x = x[keep], y = yt[keep])
  f <- if (is.finite(ub) && lb > 0) {
    frm(y | trunc(lb = 1, ub = 15) ~ x, family = discrete_weibull(),
        data = transform(d, y = y))
  } else NULL
  d$lb <- lb; d$ub <- ub
  f <- tryCatch(
    if (lb > 0 && is.finite(ub)) frm(y | trunc(lb = lb, ub = ub) ~ x,
                                     family = discrete_weibull(), data = d)
    else if (lb > 0) frm(y | trunc(lb = lb) ~ x, family = discrete_weibull(),
                         data = d)
    else frm(y | trunc(ub = ub) ~ x, family = discrete_weibull(), data = d),
    error = function(e) e)
  if (inherits(f, "error")) { cat(lb, ub, "ERROR", conditionMessage(f), "\n"); next }
  mu <- dp(f, "mu"); k <- dp(f, "shape")
  S <- function(v) ifelse(v <= 0, 1, mu^(pmax(v, 0)^k))
  Fc <- function(v) 1 - S(v + 1)
  lo <- if (lb > 0) Fc(lb - 1) else 0
  hi <- if (is.finite(ub)) Fc(ub) else 1
  ll <- log(S(d$y) - S(d$y + 1)) - log(hi - lo)
  a <- as.numeric(logLik(f))
  # candidate readings of the upper bound
  hi2 <- if (is.finite(ub)) Fc(ub - 1) else 1
  hi3 <- if (is.finite(ub)) Fc(ub + 1) else 1
  cat(sprintf("trunc(%g, %g): logLik %.10f hand %.10f diff %.3e | with F(ub - 1) %.10f | F(ub + 1) %.10f\n",
              lb, ub, a, sum(ll), a - sum(ll),
              sum(log(S(d$y) - S(d$y + 1)) - log(hi2 - lo)),
              sum(log(S(d$y) - S(d$y + 1)) - log(hi3 - lo))))
  # poisson control with the same window
  lam0 <- exp(1 + 0.2 * x)
  yp <- rpois(n, lam0)
  kp <- yp >= lb & yp <= ub
  dpp <- data.frame(x = x[kp], y = yp[kp], lb = lb, ub = ub)
  fp <- if (lb > 0 && is.finite(ub)) frm(y | trunc(lb = lb, ub = ub) ~ x,
                                         family = poisson(), data = dpp)
        else if (lb > 0) frm(y | trunc(lb = lb) ~ x, family = poisson(), data = dpp)
        else frm(y | trunc(ub = ub) ~ x, family = poisson(), data = dpp)
  lam <- dp(fp, "mu")
  lo <- if (lb > 0) ppois(lb - 1, lam) else 0
  hi <- if (is.finite(ub)) ppois(ub, lam) else 1
  llp <- dpois(dpp$y, lam, log = TRUE) - log(hi - lo)
  cat(sprintf("   poisson control: logLik %.10f hand %.10f diff %.3e\n",
              as.numeric(logLik(fp)), sum(llp), as.numeric(logLik(fp)) - sum(llp)))
}
fam <- discrete_weibull()
cat("lcdf at q = 15, mu 0.7, shape 1.1:",
    fam$lcdf(15, list(mu = 0.7, shape = 1.1, .eta_mu = qlogis(0.7)), list()),
    " hand F(15):", 1 - 0.7^(16^1.1), "\n")
f <- frm(y | trunc(lb = 2, ub = 9) ~ x, family = discrete_weibull(),
         data = data.frame(x = x, y = yt)[yt >= 2 & yt <= 9, ])
cat("trunc(2, 9) fit: mu by link", head(dp(f, "mu"), 3), "; by type = 'response'",
    head(dp_resp(f, "mu"), 3), "\n")
cat("shape by link", head(dp(f, "shape"), 2), "; response", head(dp_resp(f, "shape"), 2), "\n")
f0 <- frm(y ~ x, family = discrete_weibull(), data = data.frame(x = x, y = yt))
cat("untruncated fit: mu by link", head(dp(f0, "mu"), 3), "; response",
    head(dp_resp(f0, "mu"), 3), "\n")
