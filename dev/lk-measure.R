# Saturation measurement: for each candidate link, where does the PLAIN
# round trip log(mu) / log(1 - mu) first lose 1e-8 relative accuracy
# against a trusted reference, and where does the robust field?
options(digits = 17)
lsa <- function(a, b) RTMB::logspace_add(a, b)
sp  <- function(x) lsa(0 * x, x)                 # log1p_exp
lil <- function(x) -lsa(0 * x, -x)               # log(plogis(x))
l1l <- function(x) -lsa(0 * x, x)                # log(1 - plogis(x))

cand <- list(
  probit = list(
    inv = function(e) stats::pnorm(e),
    lo  = function(e) log(stats::pnorm(e)) - log(stats::pnorm(-e)),
    rmu = function(e) stats::pnorm(e, log.p = TRUE),
    r1m = function(e) stats::pnorm(-e, log.p = TRUE)),
  probit_approx = list(
    inv = function(e) 1 / (1 + exp(-(0.07056 * e^3 + 1.5976 * e))),
    lo  = function(e) 0.07056 * e^3 + 1.5976 * e,
    rmu = function(e) lil(0.07056 * e^3 + 1.5976 * e),
    r1m = function(e) l1l(0.07056 * e^3 + 1.5976 * e)),
  cauchit = list(
    inv = function(e) 0.5 + atan(e) / pi,
    lo  = function(e) log(0.5 + atan(e) / pi) - log(0.5 - atan(e) / pi),
    rmu = function(e) stats::pcauchy(e, log.p = TRUE),
    r1m = function(e) stats::pcauchy(e, lower.tail = FALSE, log.p = TRUE)),
  softit = list(
    inv = function(e) { y <- sp(e); y / (1 + y) },
    lo  = function(e) { y <- sp(e); log(y) - log1p(y) },
    rmu = function(e) { y <- sp(e); log(y) - log1p(y) - lsa(0*e, -(log(y) - log1p(y))) + lsa(0*e, -(log(y)-log1p(y))) - lsa(0*e, -(log(y)-log1p(y))) },
    r1m = function(e) -log1p(sp(e))),
  cloglog = list(
    inv = function(e) 1 - exp(-exp(e)),
    lo  = function(e) { t <- exp(e); log(-expm1(-t)) + t },
    rmu = function(e) log(-expm1(-exp(e))),
    r1m = function(e) -exp(e))
)
# reference log(mu) for softit done properly
cand$softit$rmu <- function(e) { y <- sp(e); log(y) - log1p(y) }

grid <- c(seq(0.5, 12, by = 0.1), seq(12.5, 60, by = 0.5),
          10^seq(2, 17, by = 0.25))
firstbad <- function(f, ref, xs, tol = 1e-8) {
  for (x in xs) {
    a <- suppressWarnings(f(x)); b <- ref(x)
    if (!is.finite(a) || abs(a - b) / max(1e-300, abs(b)) > tol) return(x)
  }
  NA_real_
}
cat(sprintf("%-14s %-22s %14s %14s\n", "link", "quantity", "plain", "robust"))
for (nm in names(cand)) {
  k <- cand[[nm]]
  # upper tail: log(1 - mu) at +eta ; lower tail: log(mu) at -eta
  p1m <- function(e) log1p(-k$inv(e))
  r1m <- k$r1m
  pmu <- function(e) log(k$inv(-e))
  rmu <- function(e) k$rmu(-e)
  cat(sprintf("%-14s %-22s %14.6g %14.6g\n", nm, "log(1-mu), eta>0",
              firstbad(p1m, r1m, grid), firstbad(r1m, r1m, grid)))
  cat(sprintf("%-14s %-22s %14.6g %14.6g\n", nm, "log(mu), eta<0",
              firstbad(pmu, rmu, grid), firstbad(rmu, rmu, grid)))
}
