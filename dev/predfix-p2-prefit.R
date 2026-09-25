# Punch round 2, B1: what does a pre-fit that should not seed the
# reported fit look like? The pre-fit is the plain fit of the model with
# the planned columns standardized, so fit exactly that (x / sd(x),
# autoscale = FALSE) and record the candidate criteria: an error, the
# code, max|grad|, the outer Hessian's eigenvalues (optimHess on the
# objective), max |par|, and the warnings the fit raised.
#   PREDFIX_ARM=lane Rscript dev/predfix-p2-prefit.R > dev/predfix-log/p2-prefit.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
probe <- function(label, fo, fam, d) {
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(fo, family = fam, data = d,
                 control = frmtmb_control(autoscale = FALSE)),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  if (inherits(f, "error")) {
    cat(sprintf("%-28s ERROR %s\n", label, substr(conditionMessage(f), 1, 60)))
    return(invisible())
  }
  H <- tryCatch(stats::optimHess(f$opt$par, f$obj$fn, f$obj$gr),
                error = function(e) NULL)
  ev <- if (is.null(H) || any(!is.finite(H))) c(NA, NA) else
    range(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  cb <- names(f$opt$par) %in% c("beta", "betad")
  Hb <- if (is.null(H) || !any(cb)) NULL else H[cb, cb, drop = FALSE]
  evb <- if (is.null(Hb) || any(!is.finite(Hb))) c(NA, NA) else
    range(eigen((Hb + t(Hb)) / 2, symmetric = TRUE,
                only.values = TRUE)$values)
  gb <- if (any(cb)) max(abs(f$obj$gr(f$opt$par)[cb])) else NA
  cat(sprintf("%-26s code %d eig [%9.2e, %9.2e] coef eig [%9.2e, %9.2e] coef|grad| %8.2e max|coef| %8.3g warn %d\n",
              label, f$opt$convergence, ev[1], ev[2], evb[1], evb[2], gb,
              if (any(cb)) max(abs(f$opt$par[cb])) else NA, length(w)))
}
cat("== reviewer design A: separated bernoulli, x sd 1e-4, standardized\n")
for (seed in 511:520) {
  set.seed(seed)
  n <- 240
  d <- data.frame(x = rnorm(n) * 1e-4, z = rnorm(n))
  d$yb <- as.integer(d$z > 0)
  d$x <- d$x / sd(d$x)
  probe(paste("A seed", seed), bf(yb ~ z + x), bernoulli(), d)
}
cat("== reviewer design B: separated z, slope on x sd 0.03, standardized\n")
for (seed in 531:540) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 12))
  d2 <- data.frame(g, x = rnorm(240) * 0.03, z = rnorm(240))
  d2$yb <- as.integer(d2$z > 0)
  d2$x <- d2$x / sd(d2$x)
  probe(paste("B seed", seed), bf(yb ~ z + x + (1 + x | g)), bernoulli(), d2)
}
cat("== the battery's random slope with no slope variance, standardized\n")
set.seed(20260923)
nb <- 200
xb <- rnorm(nb)
zb <- runif(nb)
gb <- factor(rep(1:10, length.out = nb))
fb <- factor(sample(c("a", "b", "c"), nb, TRUE))
ub <- rnorm(10, 0, 0.4)[gb]
db <- data.frame(x = xb, g = gb, gau = 0.3 + 0.4 * xb + ub + rnorm(nb))
probe("battery gau slope", bf(gau ~ x + (1 + x | g)), gaussian(), db)
cat("== healthy engaged designs, standardized (what a good pre-fit is)\n")
for (seed in 1:5) {
  set.seed(seed)
  n <- 250
  x <- rnorm(n)
  d <- data.frame(x = x, y = rpois(n, exp(0.5 + 0.4 * x)),
                  b = rbinom(n, 1, plogis(0.3 + 0.8 * x)))
  probe(paste("pois 0 + x seed", seed), bf(y ~ 0 + x), poisson(), d)
  probe(paste("bern 1 + x seed", seed), bf(b ~ x), bernoulli(), d)
  g <- factor(rep(1:20, each = 15))
  x <- rnorm(300)
  ds <- data.frame(g = g, x = x, y = 1 + 0.5 * x + rnorm(20, 0, 0.7)[g] +
                     rnorm(20, 0, 0.4)[g] * x + rnorm(300))
  probe(paste("gau slope seed", seed), bf(y ~ x + (1 + x | g)), gaussian(), ds)
  ds$y0 <- 1 + 0.5 * x + rnorm(20, 0, 0.7)[g] + rnorm(300)
  probe(paste("gau zero-var slope seed", seed), bf(y0 ~ x + (1 + x | g)),
        gaussian(), ds)
}
