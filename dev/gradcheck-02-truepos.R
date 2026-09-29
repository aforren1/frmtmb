# The fits the gradient check EXISTS for. Every candidate criterion must
# still fire on these, or the replacement is less diagnostic than the
# default it replaces.
#
#   Rscript dev/gradcheck-02-truepos.R base|lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "base"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

rows <- list()
add <- function(r) {
  rows[[length(rows) + 1L]] <<- r
  gc_print(r)
}
hr <- function(s) cat("\n== ", s, "\n", sep = "")

## ------------------------------------------------------------------
hr("T1. stopped at the iteration cap (gaussian GLM, cap 2)")
set.seed(201)
n <- 500
dd <- data.frame(x = rnorm(n), z = rnorm(n))
dd$y <- rnorm(n, 3 + 2 * dd$x - 1.5 * dd$z, 1)
ref <- frm(bf(y ~ x + z), family = gaussian(), data = dd)
for (cap in c(1L, 2L, 3L, 5L, 10L)) {
  cp <- gc_catch(frm(bf(y ~ x + z), family = gaussian(), data = dd,
                     control = frmtmb_control(
                       restarts = 0,
                       optCtrl = list(iter.max = cap, eval.max = cap * 3))))
  f <- cp$value
  add(gc_row(paste0("gaussian, iter.max = ", cap), 201, f, cp$warnings))
  cat("  logLik shortfall against the converged fit:",
      format(as.numeric(logLik(ref)) - as.numeric(logLik(f)),
             digits = 6), "\n")
}

## ------------------------------------------------------------------
hr("T2. started far away with a tiny iteration cap")
for (cap in c(1L, 2L, 4L, 8L)) {
  cp <- gc_catch(frm(bf(y ~ x + z), family = gaussian(), data = dd,
                     start = list(beta = c(50, -40, 30)),
                     control = frmtmb_control(
                       restarts = 0,
                       optCtrl = list(iter.max = cap, eval.max = cap * 3))))
  f <- cp$value
  add(gc_row(paste0("far start, iter.max = ", cap), 201, f, cp$warnings))
  cat("  logLik shortfall:",
      format(as.numeric(logLik(ref)) - as.numeric(logLik(f)),
             digits = 6), "\n")
}

## ------------------------------------------------------------------
hr("T3. a poisson GLMM stopped at the cap")
set.seed(202)
n <- 600
gd <- data.frame(x = rnorm(n), g = factor(rep(1:30, 20)))
gd$y <- rpois(n, exp(0.5 + 0.6 * gd$x + rnorm(30, 0, 0.7)[gd$g]))
refg <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = gd)
for (cap in c(1L, 2L, 4L, 8L)) {
  cp <- gc_catch(frm(bf(y ~ x + (1 | g)), family = poisson(), data = gd,
                     control = frmtmb_control(
                       restarts = 0,
                       optCtrl = list(iter.max = cap, eval.max = cap * 3))))
  f <- cp$value
  add(gc_row(paste0("poisson GLMM, iter.max = ", cap), 202, f, cp$warnings))
  cat("  logLik shortfall:",
      format(as.numeric(logLik(refg)) - as.numeric(logLik(f)),
             digits = 6), "\n")
}

## ------------------------------------------------------------------
hr("T4. a flat direction: a bump whose centre starts off the data")
# The nonlinear peak fit from diagnose_flat()'s own docstring: the
# likelihood does not depend on lamp, pk or lsig at all at the start.
set.seed(203)
nw <- 45
pd <- data.frame(w = seq_len(nw))
pd$y <- rnorm(nw, 2, 0.5)
cp <- gc_catch(frm(
  bf(y ~ exp(lamp) * exp(-0.5 * ((w - pk) / exp(lsig))^2),
     nlf(lamp ~ 1), nlf(pk ~ 1), nlf(lsig ~ 1), nl = TRUE),
  family = gaussian(), data = pd))
f <- cp$value
if (inherits(f, "gc_error")) {
  cat("  ERROR:", f, "\n")
} else {
  add(gc_row("flat nonlinear peak", 203, f, cp$warnings))
  cat("  warnings:\n"); for (w in cp$warnings) cat("   -", w, "\n")
}

## ------------------------------------------------------------------
hr("T5. a genuinely unconverged ordinal fit (cap 2, n = 6000)")
set.seed(204)
n <- 6000
od <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
eo <- 0.8 * od$x1 - 0.5 * od$x2
od$yo <- cut(eo + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
             labels = FALSE)
refo <- frm(bf(yo ~ x1 + x2), family = cumulative(), data = od)
for (cap in c(2L, 4L, 8L)) {
  cp <- gc_catch(frm(bf(yo ~ x1 + x2), family = cumulative(), data = od,
                     control = frmtmb_control(
                       restarts = 0,
                       optCtrl = list(iter.max = cap, eval.max = cap * 3))))
  f <- cp$value
  add(gc_row(paste0("cumulative n=6000, cap ", cap), 204, f, cp$warnings))
  cat("  logLik shortfall:",
      format(as.numeric(logLik(refo)) - as.numeric(logLik(f)),
             digits = 6), "\n")
}

## ------------------------------------------------------------------
hr("T6. a bound that is NOT active but the fit is unconverged")
cp <- gc_catch(frm(bf(y ~ x + z), family = gaussian(), data = dd,
                   prior = set_prior("", class = "b", ub = 100),
                   start = list(beta = c(30, -20, 20)),
                   control = frmtmb_control(
                     restarts = 0,
                     optCtrl = list(iter.max = 2, eval.max = 6))))
f <- cp$value
add(gc_row("inactive bound, cap 2", 201, f, cp$warnings))

## ------------------------------------------------------------------
hr("T7. cost of the Hessian, which the new criterion needs")
timeit <- function(lbl, f) {
  p <- f$opt$par
  t1 <- proc.time()[["elapsed"]]
  reps <- 0L
  repeat {
    H <- stats::optimHess(p, function(q) f$obj$fn(q),
                          function(q) drop(f$obj$gr(q)))
    reps <- reps + 1L
    if (proc.time()[["elapsed"]] - t1 > 1.5) break
  }
  el <- (proc.time()[["elapsed"]] - t1) / reps
  t2 <- proc.time()[["elapsed"]]
  g <- 0
  gr <- 0L
  repeat {
    g <- f$obj$gr(p)
    gr <- gr + 1L
    if (proc.time()[["elapsed"]] - t2 > 1.0) break
  }
  elg <- (proc.time()[["elapsed"]] - t2) / gr
  cat(sprintf("%-28s np %3d  optimHess %8.4f s  one gr %8.5f s  ratio %6.1f\n",
              lbl, length(p), el, elg, el / elg))
}
timeit("gaussian GLM n=500", ref)
timeit("poisson GLMM n=600 q=30", refg)
timeit("cumulative n=6000", refo)

set.seed(205)
nb <- 4000
bigd <- data.frame(g = factor(rep(1:80, each = 50)))
bigd$fx <- factor(rep(1:25, length.out = nb))
bigd$x <- rnorm(nb)
bigd$y <- rnorm(nb, 1 + 0.4 * bigd$x + rnorm(80, 0, 0.6)[bigd$g] +
                 rnorm(25, 0, 0.5)[bigd$fx], 1)
bigf <- frm(bf(y ~ x + fx + (1 | g)), family = gaussian(), data = bigd)
timeit("gaussian LMM 27 fixed pars", bigf)

saveRDS(rows, file.path("dev", paste0("gradcheck-02-", which_lib, ".rds")))
cat("\nrows:", length(rows), "  warned:",
    sum(vapply(rows, function(r) isTRUE(r$warned), TRUE)), "\n")
