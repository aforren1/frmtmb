# Reviewer probe: the discrete Weibull mean against a reference that sums
# 1e7 terms before the incomplete-gamma tail, so the Euler-Maclaurin
# cut at m = 1000 that dw_moments() uses is tested; and GEV fits at
# response scales 1e-3 and 1e3 (the predictor unit the lane added).
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
dwm <- frmtmb:::dw_moments
ref <- function(mu, k, M) {
  L <- log(mu)
  y <- seq_len(M)
  s <- sum(exp(L * y^k))
  c0 <- M + 1
  f0 <- exp(L * c0^k)
  a <- 1 / k
  x0 <- -L * c0^k
  ti <- exp(lgamma(a) + pgamma(x0, a, lower.tail = FALSE, log.p = TRUE) -
              log(k) - a * log(-L))
  d0 <- f0 * L * k * c0^(k - 1)
  s + ti + f0 / 2 - d0 / 12
}
cat("=== A. dw mean, m = 1000 against M = 1e5 and 1e7\n")
g <- expand.grid(mu = c(0.3, 0.9, 0.99, 0.999), shape = c(0.1, 0.2, 0.35, 0.5))
for (i in seq_len(nrow(g))) {
  m <- dwm(g$mu[i], g$shape[i])$mean
  r5 <- ref(g$mu[i], g$shape[i], 1e5)
  r7 <- ref(g$mu[i], g$shape[i], 1e7)
  cat(sprintf("mu %.3f shape %.2f: dw_moments %.12g  M=1e5 %.12g  M=1e7 %.12g  rel(m, 1e7) %.2e  rel(1e5, 1e7) %.2e\n",
              g$mu[i], g$shape[i], m, r5, r7, abs(m - r7) / r7,
              abs(r5 - r7) / r7))
}

cat("\n=== B. GEV fits at response scales 1e-3, 1, 1e3 (seed 1201)\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
set.seed(1201)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(15, 20))
b <- rnorm(15, 0, 0.4)[d$g]
y0 <- rg(n, 1 + 0.5 * d$x + b, exp(-0.2 + 0.2 * d$z), 0.15)
for (sc in c(1e-3, 1, 1e3)) {
  d$y <- y0 * sc
  w <- character(0)
  f <- withCallingHandlers(frm(bf(y ~ x + (1 | g), sigma ~ z),
                               family = gen_extreme_value(), data = d),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  xi <- summary(f)$spec_pars["xi", 1]
  cat(sprintf("scale %g: conv %d logLik + n log(scale) %.8f xi %.8f sd(g)/scale %.8f warnings %d\n",
              sc, f$opt$convergence, as.numeric(logLik(f)) + n * log(sc), xi,
              VarCorr(f)[["g"]][["sd"]][1, "Estimate"] / sc, length(w)))
}
