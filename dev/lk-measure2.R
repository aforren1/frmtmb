# Saturation, measured properly.
#   plain  : mu <- linkinv(eta); log(mu), log1p(-mu)
#   robust : lo <- logit_eta(eta); log_inv_logit(lo), log1m_inv_logit(lo)
#            or lm <- log_eta(eta) directly
#   ref    : an independent accurate closed form
# Reported: the first |eta| on the grid at which the path's relative
# error exceeds 1e-8 or it stops being finite.
options(digits = 17)
lsa <- function(a, b) RTMB::logspace_add(a, b)
sp  <- function(x) lsa(0 * x, x)
lil <- function(x) -lsa(0 * x, -x)
l1l <- function(x) -lsa(0 * x, x)

L <- list(
  probit = list(
    inv = function(e) stats::pnorm(e),
    lo  = function(e) log(stats::pnorm(e)) - log(stats::pnorm(-e)),
    ref_mu  = function(e) stats::pnorm(e, log.p = TRUE),
    ref_1m  = function(e) stats::pnorm(-e, log.p = TRUE)),
  probit_approx = list(
    inv = function(e) 1 / (1 + exp(-(0.07056 * e^3 + 1.5976 * e))),
    lo  = function(e) 0.07056 * e^3 + 1.5976 * e,
    ref_mu = function(e) stats::plogis(0.07056 * e^3 + 1.5976 * e, log.p = TRUE),
    ref_1m = function(e) stats::plogis(0.07056 * e^3 + 1.5976 * e,
                                       lower.tail = FALSE, log.p = TRUE)),
  cauchit = list(
    inv = function(e) 0.5 + atan(e) / pi,
    lo  = function(e) log(0.5 + atan(e) / pi) - log(0.5 - atan(e) / pi),
    ref_mu = function(e) stats::pcauchy(e, log.p = TRUE),
    ref_1m = function(e) stats::pcauchy(e, lower.tail = FALSE, log.p = TRUE)),
  softit = list(
    inv = function(e) { y <- sp(e); y / (1 + y) },
    lo  = function(e) { y <- sp(e); log(y) - log1p(y) },
    # reference from the asymptote: y = log1p(exp(e)) with log1p/expm1
    ref_mu = function(e) { y <- log1p(exp(e)); if (e < -30) y <- exp(e) - exp(2*e)/2 + exp(3*e)/3
                           log(y) - log1p(y) },
    ref_1m = function(e) { y <- if (e > 30) e + log1p(exp(-e)) else log1p(exp(e))
                           -log1p(y) }),
  cloglog = list(
    inv = function(e) 1 - exp(-exp(e)),
    lo  = function(e) { t <- exp(e); log(-expm1(-t)) + t },
    ref_mu = function(e) log(-expm1(-exp(e))),
    ref_1m = function(e) -exp(e))
)
grid <- c(seq(0.5, 12, by = 0.1), seq(12.2, 80, by = 0.2),
          10^seq(2, 17, by = 0.25))
firstbad <- function(f, ref, tol = 1e-8) {
  for (x in grid) {
    a <- suppressWarnings(f(x)); b <- ref(x)
    if (!is.finite(a) || !is.finite(b)) return(if (is.finite(b)) x else NA_real_)
    if (abs(a - b) / max(1e-300, abs(b)) > tol) return(x)
  }
  NA_real_
}
cat(sprintf("%-14s %-18s %13s %13s\n", "link", "quantity", "plain", "robust"))
for (nm in names(L)) {
  k <- L[[nm]]
  cat(sprintf("%-14s %-18s %13.5g %13.5g\n", nm, "log(1-mu) e>0",
    firstbad(function(e) log1p(-k$inv(e)), k$ref_1m),
    firstbad(function(e) l1l(k$lo(e)), k$ref_1m)))
  cat(sprintf("%-14s %-18s %13.5g %13.5g\n", nm, "log(mu) e<0",
    firstbad(function(e) log(k$inv(-e)), function(e) k$ref_mu(-e)),
    firstbad(function(e) lil(k$lo(-e)), function(e) k$ref_mu(-e))))
}
cat("\n-- positive-mean links, log(mu) --\n")
P <- list(
  softplus = list(inv = function(e) sp(e), rob = NULL,
    ref = function(e) if (e < -30) e + log1p(-exp(e)/2 + exp(2*e)/3) else log(log1p(exp(e)))),
  squareplus = list(inv = function(e) (e + sqrt(e^2 + 4)) / 2,
    rob = function(e) asinh(e / 2),
    ref = function(e) if (e < -1e6) -log(-e) + log1p(-1/e^2) else log((e + sqrt(e^2+4))/2)),
  sqrt = list(inv = function(e) e^2, rob = function(e) 2 * log(abs(e)),
    ref = function(e) 2 * log(abs(e))),
  invsq = list(inv = function(e) 1 / sqrt(e), rob = function(e) -0.5 * log(e),
    ref = function(e) -0.5 * log(e))
)
gridn <- c(-10^seq(0, 17, by = 0.25), 10^seq(-17, 17, by = 0.25))
fb2 <- function(f, ref, xs, tol = 1e-8) {
  for (x in xs) {
    a <- suppressWarnings(f(x)); b <- suppressWarnings(ref(x))
    if (!is.finite(b)) next
    if (!is.finite(a) || abs(a - b) / max(1e-300, abs(b)) > tol) return(x)
  }
  NA_real_
}
for (nm in names(P)) {
  k <- P[[nm]]
  xs <- if (nm %in% c("softplus", "squareplus", "sqrt")) gridn else 10^seq(-300, 300, by = 5)
  if (nm == "sqrt") xs <- c(-10^seq(0,-160,by=-2), 10^seq(0,-160,by=-2))
  cat(sprintf("%-12s plain %13.5g   robust %13.5g\n", nm,
              fb2(function(e) log(k$inv(e)), k$ref, xs),
              if (is.null(k$rob)) NA_real_ else fb2(k$rob, k$ref, xs)))
}
