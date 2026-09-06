options(digits = 17)
lsa <- function(a, b) RTMB::logspace_add(a, b)
sp  <- function(x) lsa(0 * x, x)
lil <- function(x) -lsa(0 * x, -x)
l1l <- function(x) -lsa(0 * x, x)
cat("-- softit: mu = y/(1+y), y = softplus(eta); log-odds is log(y) --\n")
for (e in c(-3, -0.5, 0, 1, 4)) {
  y <- log1p(exp(e)); mu <- y / (1 + y)
  cat(sprintf("eta %5.1f  logit(mu)=%.16g  log(y)=%.16g\n", e,
              log(mu / (1 - mu)), log(y)))
}
grid <- c(seq(0.5, 12, by = 0.1), seq(12.2, 80, by = 0.2), 10^seq(2, 17, by = 0.25))
fb <- function(f, ref, xs = grid, tol = 1e-8) {
  for (x in xs) {
    a <- suppressWarnings(f(x)); b <- suppressWarnings(ref(x))
    if (!is.finite(b)) return(NA_real_)
    if (!is.finite(a) || abs(a - b) / max(1e-300, abs(b)) > tol) return(x)
  }
  NA_real_
}
inv <- function(e) { y <- sp(e); y / (1 + y) }
lo  <- function(e) log(sp(e))
ref1m <- function(e) -log1p(if (e > 30) e + log1p(exp(-e)) else log1p(exp(e)))
refmu <- function(e) { y <- if (e < -30) exp(e) - exp(2*e)/2 + exp(3*e)/3 else log1p(exp(e))
                       log(y) - log1p(y) }
cat(sprintf("softit log(1-mu) e>0  plain %11.5g  robust %11.5g\n",
            fb(function(e) log1p(-inv(e)), ref1m),
            fb(function(e) l1l(lo(e)), ref1m)))
cat(sprintf("softit log(mu)   e<0  plain %11.5g  robust %11.5g\n",
            fb(function(e) log(inv(-e)), function(e) refmu(-e)),
            fb(function(e) lil(lo(-e)), function(e) refmu(-e))))
cat("\n-- squareplus: exact log(mu) is asinh(eta/2) --\n")
for (e in c(-8, -1, 0, 1, 8)) {
  cat(sprintf("eta %5.1f  log(linkinv)=%.16g  asinh(e/2)=%.16g\n", e,
              log((e + sqrt(e^2 + 4)) / 2), asinh(e / 2)))
}
neg <- -10^seq(0, 12, by = 0.25)
cat(sprintf("squareplus log(mu) eta<0  plain %11.5g (vs asinh)\n",
            fb(function(e) log((e + sqrt(e^2 + 4)) / 2),
               function(e) asinh(e / 2), neg)))
cat(sprintf("  asinh vs the -log(-eta) asymptote, eta < -1e6: max rel %.3g\n",
            max(abs(asinh(neg[neg < -1e6]/2) - (-log(-neg[neg < -1e6]))) /
                abs(log(-neg[neg < -1e6])))))
cat("\n-- softplus log(mu): where does log(logspace_add(0,eta)) die --\n")
for (e in c(-700, -740, -745, -746, -750)) {
  cat(sprintf("eta %7.1f  log(sp(eta)) = %.10g\n", e, log(sp(e))))
}
