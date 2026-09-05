# Self-diagnosing version. Three models, easiest first, so that a NaN can
# be attributed: (A) scalar intercept, the case frmtmb runs today;
# (B) intercept + slope, INDEPENDENT (diag), prior as a plain vectorized
# dnorm; (C) intercept + slope, correlated (us). Each is taped at two
# anchors, the Laplace optimum and the cold template, as quad_fit() does.
suppressMessages(library(RTMB))
set.seed(11)
G <- 60L; m <- 8L
g <- rep(seq_len(G), each = m)
x <- rnorm(G * m)
Sig <- matrix(c(1, 0.35, 0.35, 0.49), 2)
B <- MASS::mvrnorm(G, c(0, 0), Sig)
y <- rbinom(G * m, 1, plogis(-0.5 + 0.8 * x + B[g, 1] + B[g, 2] * x))

models <- list(
  A = list(
    tpl = list(beta = c(0, 0), theta = 0, b = numeric(G)), dim = 1L,
    nll = function(p) {
      getAll(p)
      et <- beta[1] + beta[2] * x + b[g]
      -sum(dnorm(b, 0, exp(theta), log = TRUE)) -
        sum(dbinom(y, 1, plogis(et), log = TRUE))
    }),
  B = list(
    tpl = list(beta = c(0, 0), theta = c(0, 0), b = numeric(2 * G)), dim = 2L,
    nll = function(p) {
      getAll(p)
      bi <- b[seq(1, 2 * G, by = 2)]; bs <- b[seq(2, 2 * G, by = 2)]
      et <- beta[1] + beta[2] * x + bi[g] + bs[g] * x
      -sum(dnorm(bi, 0, exp(theta[1]), log = TRUE)) -
        sum(dnorm(bs, 0, exp(theta[2]), log = TRUE)) -
        sum(dbinom(y, 1, plogis(et), log = TRUE))
    }),
  C = list(
    tpl = list(beta = c(0, 0), theta = c(0, 0, 0), b = numeric(2 * G)), dim = 2L,
    nll = function(p) {
      getAll(p)
      sd <- exp(theta[1:2]); rho <- tanh(theta[3])
      bi <- b[seq(1, 2 * G, by = 2)]; bs <- b[seq(2, 2 * G, by = 2)]
      # the us prior written as a Cholesky whitening, so the density is a
      # vectorized sum of dnorm terms plus the log-determinant, with no
      # per-group loop and no matrix reshaping of b
      z1 <- bi / sd[1]
      z2 <- (bs - rho * sd[2] * z1) / (sd[2] * sqrt(1 - rho^2))
      lp <- sum(dnorm(z1, log = TRUE)) + sum(dnorm(z2, log = TRUE)) -
        G * (log(sd[1]) + log(sd[2]) + 0.5 * log(1 - rho^2))
      et <- beta[1] + beta[2] * x + bi[g] + bs[g] * x
      -lp - sum(dbinom(y, 1, plogis(et), log = TRUE))
    })
)

gk <- function(dim, adaptive, debug = FALSE) list(b = structure(
  list(dim = dim, adaptive = adaptive, debug = debug, method = "marginal_gk"),
  class = "GK"))

for (nm in names(models)) {
  M <- models[[nm]]
  lap <- MakeADFun(M$nll, M$tpl, random = "b", silent = TRUE)
  lo <- nlminb(lap$par, lap$fn, lap$gr)
  cat(sprintf("\n[%s] Laplace nll %.4f  par %s\n", nm, lo$objective,
              paste(round(lo$par, 3), collapse = " ")))
  anchors <- list(laplace = lap$env$parList(lo$par, lap$env$last.par.best),
                  cold = M$tpl)
  for (an in names(anchors)) for (ad in c(FALSE, TRUE)) {
    t0 <- Sys.time()
    o <- try(MakeADFun(M$nll, anchors[[an]], random = "b",
                       integrate = gk(M$dim, ad), silent = TRUE), silent = TRUE)
    if (inherits(o, "try-error")) {
      cat(sprintf("  GK dim=%d adaptive=%-5s anchor=%-7s: MakeADFun error: %s\n",
                  M$dim, ad, an, trimws(conditionMessage(attr(o, "condition")))))
      next
    }
    f0 <- try(o$fn(o$par), silent = TRUE)
    if (inherits(f0, "try-error") || !is.finite(f0)) {
      cat(sprintf("  GK dim=%d adaptive=%-5s anchor=%-7s: fn at anchor = %s (tape %.1fs)\n",
                  M$dim, ad, an, if (inherits(f0, "try-error")) "error" else f0,
                  as.numeric(Sys.time() - t0, units = "secs")))
      next
    }
    t1 <- Sys.time()
    op <- try(nlminb(o$par, o$fn, o$gr), silent = TRUE)
    if (inherits(op, "try-error")) {
      cat(sprintf("  GK dim=%d adaptive=%-5s anchor=%-7s: fn %.4f at anchor, optimizer error\n",
                  M$dim, ad, an, f0)); next
    }
    cat(sprintf("  GK dim=%d adaptive=%-5s anchor=%-7s: nll %.4f  par %s  (opt %.1fs, %d fn evals)\n",
                M$dim, ad, an, op$objective, paste(round(op$par, 3), collapse = " "),
                as.numeric(Sys.time() - t1, units = "secs"),
                op$evaluations[["function"]]))
  }
}
# one debug run on the smallest failing case, to see what the transform says
cat("\n== transform debug output, model A, adaptive=FALSE, laplace anchor ==\n")
M <- models$A
lap <- MakeADFun(M$nll, M$tpl, random = "b", silent = TRUE)
lo <- nlminb(lap$par, lap$fn, lap$gr)
o <- try(MakeADFun(M$nll, lap$env$parList(lo$par, lap$env$last.par.best),
                   random = "b", integrate = gk(1L, FALSE, debug = TRUE),
                   silent = TRUE))
if (!inherits(o, "try-error")) print(o$fn(o$par))
if (requireNamespace("GLMMadaptive", quietly = TRUE)) {
  d <- data.frame(y = y, x = x, g = factor(g))
  ga <- GLMMadaptive::mixed_model(y ~ x, ~ x | g, d, binomial(), nAGQ = 15)
  cat(sprintf("\nGLMMadaptive nAGQ=15 (model C reference): nll %.4f  beta %s\n",
              -as.numeric(logLik(ga)),
              paste(round(GLMMadaptive::fixef(ga), 3), collapse = " ")))
}
