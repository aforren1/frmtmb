## Reviewer, claim 3 (the misses) and the importance channel.
## For each row: the fit, a reference fit of the same data optimized
## hard, the true shortfall, the headroom the check measured, and
## whether the build warned. A MISS is shortfall > grad_tol with no
## warning of any kind.
## usage: Rscript gradcheck-rev-09-miss.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")
`%||%` <- function(x, y) if (is.null(x)) y else x
gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                             warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w,
       grad = any(grepl("Large maximum absolute gradient", w, fixed = TRUE)))
}
tight <- function() frmtmb_control(
  restarts = 8,
  optCtrl = list(rel.tol = 1e-14, x.tol = 1e-14, iter.max = 5000,
                 eval.max = 5000))
loose <- function(rt) frmtmb_control(
  restarts = 0,
  optCtrl = list(rel.tol = rt, x.tol = rt, iter.max = 2000,
                 eval.max = 2000))

show <- function(tag, r, ref) {
  f <- r$fit
  if (inherits(f, "condition")) {
    cat(sprintf("%-34s ERROR %s\n", tag, conditionMessage(f))); return(NULL)
  }
  d <- tryCatch(diagnose(f, quiet = TRUE), error = function(e) NULL)
  tol <- f$control$grad_tol
  sh <- if (is.null(ref) || inherits(ref, "condition")) NA_real_ else
    as.numeric(logLik(ref)) - as.numeric(logLik(f))
  anyw <- length(r$w) > 0
  miss <- isTRUE(sh > tol) && !anyw
  cat(sprintf(paste0("%-34s conv=%d gmax=%-10s head=%-10s shortfall=%-11s",
                     " gradwarn=%-5s anywarn=%-5s MISS=%s\n"),
              tag, f$opt$convergence,
              format(d$max_grad %||% NA, digits = 4),
              format(d$grad_headroom %||% NA, digits = 4),
              format(sh, digits = 6), r$grad, anyw, miss))
  if (!is.na(sh) && !is.null(d) && is.finite(d$grad_headroom %||% NA) &&
      (d$grad_headroom %||% 0) > 0) {
    cat(sprintf("%-34s shortfall / headroom = %s\n", "",
                format(sh / d$grad_headroom, digits = 4)))
  }
  for (w in r$w) cat("      WARN:", substr(w, 1, 110), "\n")
  invisible(d)
}

## ---- M1: the trip-wire's gap the worker recorded (defect 1), on BOTH
## builds: a column with sd 1e-7 and autoscale = FALSE.
set.seed(501)
n <- 400
d1 <- data.frame(xs = rnorm(n) * 1e-7)
d1$y <- rnorm(n, 1 + 2e7 * d1$xs, 1)
noas <- function(extra = list()) {
  do.call(frmtmb_control, c(list(autoscale = FALSE), extra))
}
r1 <- gw(frm(bf(y ~ xs), family = gaussian(), data = d1,
             control = noas()))
ref1 <- frm(bf(y ~ xs), family = gaussian(), data = d1)   # autoscale on
show("M1 sd 1e-7 col, autoscale=FALSE", r1, ref1)
cat("\n")

## ---- M2: the poisson rel.tol sweep, which holds the 2.256e-4 miss.
set.seed(404)
n <- 2000
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rpois(n, exp(0.4 + 0.7 * d2$x - 0.4 * d2$z))
ref2 <- frm(bf(y ~ x + z), family = poisson(), data = d2, control = tight())
for (rt in c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)) {
  r <- gw(frm(bf(y ~ x + z), family = poisson(), data = d2,
              control = loose(rt)))
  show(sprintf("M2 poisson rel.tol=%.0e", rt), r, ref2)
}
cat("\n")

## ---- M3: candidates for a LARGE shortfall the Newton decrement
## understates, because the Hessian is nearly singular and the valley
## is curved rather than quadratic.
## (a) a collinear ridge, at three collinearity levels
for (eps in c(1e-5, 1e-6, 1e-7, 1e-8)) {
  set.seed(502)
  n <- 1000
  x1 <- rnorm(n)
  dd <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, eps))
  dd$y <- rnorm(n, 1 + 2 * dd$x1, 1)
  r <- gw(frm(bf(y ~ x1 + x2), family = gaussian(), data = dd,
              control = frmtmb_control(restarts = 0, autoscale = FALSE)))
  ref <- gw(frm(bf(y ~ x1 + x2), family = gaussian(), data = dd,
                control = do.call(frmtmb_control,
                                  c(list(autoscale = FALSE),
                                    list(restarts = 8,
                                         optCtrl = list(rel.tol = 1e-14,
                                                        x.tol = 1e-14,
                                                        iter.max = 5000,
                                                        eval.max = 5000))))))
  show(sprintf("M3a collinear eps=%.0e", eps), r, ref$fit)
}
cat("\n")

## (b) near separation: the MLE is finite but far out, and the curvature
## there is tiny
for (rt in c(1e-2, 1e-4, 1e-6)) {
  set.seed(503)
  n <- 2000
  dd <- data.frame(x = rnorm(n))
  dd$y <- as.integer(dd$x > 0)
  dd$y[which.max(dd$x)] <- 0L
  r <- gw(frm(bf(y ~ x), family = bernoulli(), data = dd,
              control = loose(rt)))
  ref <- gw(frm(bf(y ~ x), family = bernoulli(), data = dd,
                control = tight()))
  show(sprintf("M3b near-separation rel.tol=%.0e", rt), r, ref$fit)
}
cat("\n")

## (c) a nonlinear decay stopped on its plateau, from a start far out
for (b0 in c(5, 12)) {
  set.seed(504)
  n <- 600
  dd <- data.frame(x = runif(n, 0, 5))
  dd$y <- rnorm(n, 3 * exp(-0.6 * dd$x), 0.2)
  r <- gw(frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) + gaussian(),
              data = dd, start = list(beta = c(0.2, b0)),
              control = loose(1e-8)))
  ref <- gw(frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) + gaussian(),
                data = dd, control = tight()))
  show(sprintf("M3c nl decay, b start %g", b0), r, ref$fit)
}
cat("\n")

## (d) a collapsing mixture component, stopped early
for (rt in c(1e-2, 1e-4)) {
  set.seed(505)
  n <- 1500
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x, 1)
  r <- gw(frm(bf(y ~ x), family = mixture(gaussian(), gaussian()),
              data = dd, control = loose(rt)))
  ref <- gw(frm(bf(y ~ x), family = mixture(gaussian(), gaussian()),
                data = dd, control = tight()))
  show(sprintf("M3d mixture rel.tol=%.0e", rt), r, ref$fit)
}
cat("\n")

## ---- M4: the importance channel. check_convergence() deliberately
## says nothing about the gradient of a Monte Carlo objective, but
## diagnose() now runs the same verdict on it.
set.seed(506)
ng <- 60
d4 <- data.frame(g = factor(rep(seq_len(ng), 25)))
d4$x <- rnorm(nrow(d4))
re <- rnorm(ng, 0, 0.6)
d4$y <- rpois(nrow(d4), exp(0.3 + 0.5 * d4$x + re[d4$g]))
for (rt in c(1e-3, 1e-2, 1e-1)) {
  r <- gw(frm(bf(y ~ x + (1 | g)), family = poisson(), data = d4,
              importance = 500L, control = loose(rt)))
  if (inherits(r$fit, "condition")) {
    cat("M4 rel.tol", rt, "ERROR", conditionMessage(r$fit), "\n"); next
  }
  d <- diagnose(r$fit, quiet = TRUE)
  out <- capture.output(diagnose(r$fit))
  cat(sprintf(paste0("M4 importance rel.tol=%.0e  max_grad=%-10s",
                     " grad_proj=%-10s headroom=%-10s clean=%-5s nwarn=%d\n"),
              rt, format(d$max_grad, digits = 4),
              format(d$grad_proj %||% NA, digits = 4),
              format(d$grad_headroom %||% NA, digits = 4),
              any(grepl("No convergence problems", out)), length(r$w)))
  for (w in r$w) cat("      WARN:", substr(w, 1, 110), "\n")
}
cat("DONE miss\n")
