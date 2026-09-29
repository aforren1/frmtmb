# False alarms of "Large maximum absolute gradient", constructed on the
# reference build (0.64.0) and scored by every candidate criterion.
#
#   Rscript dev/gradcheck-01-construct.R base   # rellib-r3
#   Rscript dev/gradcheck-01-construct.R lane   # wt-gradcheck-lib
#
# Each block records its own evidence that the fit is CORRECT, so the
# warning it raises is a false alarm and not a finding.

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
hr("A. an active upper bound on a fixed effect (arcov's 16.7)")
# A coefficient whose unconstrained optimum is far above the bound: the
# constrained optimum sits ON the bound, where the gradient is nonzero
# BY CONSTRUCTION. Evidence that the fit is correct: the profile is
# monotone up to the bound, so the constrained maximum is the boundary
# point, and re-optimizing from it does not move.
set.seed(101)
n <- 300
dd <- data.frame(x = rnorm(n))
dd$y <- rnorm(n, 1 + 2 * dd$x, 1)
cap <- gc_catch(frm(bf(y ~ x), family = gaussian(), data = dd,
                    prior = set_prior("", class = "b", ub = 0.1)))
fb <- cap$value
add(gc_row("bound ub=0.1 on b (gaussian)", 101, fb, cap$warnings))
cat("  estimate of x :", format(unname(coef(fb)["x"]), digits = 10), "\n")
cat("  the bound     : 0.1\n")
# the boundary IS the constrained optimum: any feasible move is worse
prof <- vapply(c(0.1, 0.099, 0.09, 0.05, 0), function(v) {
  p <- fb$opt$par
  p[["beta"]] <- p[["beta"]]  # names only
  q <- fb$opt$par
  q[which(frmtmb:::outer_par_names(fb) == "x")] <- v
  as.numeric(fb$obj$fn(q))
}, 0)
cat("  nll at x =", paste(c(0.1, 0.099, 0.09, 0.05, 0), collapse = ", "),
    ":\n    ", paste(format(prof, digits = 12), collapse = "  "), "\n")
cat("  monotone increasing away from the bound:",
    all(diff(prof) > 0), "\n")
# a restart from the optimum does not move it
re <- gc_catch(frm(bf(y ~ x), family = gaussian(), data = dd,
                   prior = set_prior("", class = "b", ub = 0.1),
                   start = list(beta = unname(fb$opt$par[1:2])),
                   control = frmtmb_control(restarts = 5)))
cat("  refit from the optimum, 5 restarts: objective delta",
    format(re$value$opt$objective - fb$opt$objective, digits = 4),
    " max|par delta| ",
    format(max(abs(re$value$opt$par - fb$opt$par)), digits = 4), "\n")

## ------------------------------------------------------------------
hr("B. an active bound on an AR coefficient (arcov's 32.9)")
set.seed(102)
ng <- 40
nt <- 15
gd <- expand.grid(t = seq_len(nt), g = factor(seq_len(ng)))
gd <- gd[order(gd$g, gd$t), ]
gd$y <- as.numeric(unlist(lapply(seq_len(ng), function(i) {
  as.numeric(arima.sim(list(ar = 0.7), nt, sd = 1)) + 1
})))
capb <- gc_catch(frm(bf(y ~ 1 + ar(time = t, gr = g, cov = TRUE)),
                     family = gaussian(), data = gd,
                     prior = set_prior("", class = "ar", ub = 0.1)))
fa <- capb$value
if (inherits(fa, "gc_error")) {
  cat("  ERROR:", fa, "\n")
} else {
  add(gc_row("bound ub=0.1 on ar (cov=TRUE)", 102, fa, capb$warnings))
  cat("  ar estimate:",
      format(fa$opt$par[frmtmb:::outer_par_names(fa) == "ar[1]"],
             digits = 10), "\n")
}

## ------------------------------------------------------------------
hr("C. a cumulative ordinal fit against MASS::polr")
# Evidence: the log-likelihood and every coefficient agree with polr,
# an independent implementation, to the digit.
set.seed(103)
n <- 6000
dd <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
eta <- 0.8 * dd$x1 - 0.5 * dd$x2
dd$yo <- cut(eta + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
             labels = FALSE)
dd$yof <- factor(dd$yo, ordered = TRUE)
capc <- gc_catch(frm(bf(yo ~ x1 + x2), family = cumulative(), data = dd))
fc <- capc$value
add(gc_row("cumulative, n = 6000", 103, fc, capc$warnings))
pl <- MASS::polr(yof ~ x1 + x2, data = dd, method = "logistic")
cat("  logLik frmtmb:", format(as.numeric(logLik(fc)), digits = 12), "\n")
cat("  logLik polr  :", format(as.numeric(logLik(pl)), digits = 12), "\n")
cat("  relative gap :",
    format(abs(as.numeric(logLik(fc)) - as.numeric(logLik(pl))) /
             abs(as.numeric(logLik(pl))), digits = 4), "\n")
cat("  coef frmtmb  :", format(coef(fc)[c("x1", "x2")], digits = 10), "\n")
cat("  coef polr    :", format(coef(pl), digits = 10), "\n")
cat("  max relative coefficient gap:",
    format(max(abs(coef(fc)[c("x1", "x2")] - coef(pl)) /
                 abs(coef(pl))), digits = 4), "\n")

## ------------------------------------------------------------------
hr("D. the same cumulative fit, re-optimized from its own optimum")
# A point the optimizer will not leave is a stationary point of the
# problem it was given, whatever the absolute gradient reads.
capd <- gc_catch(frm(bf(yo ~ x1 + x2), family = cumulative(), data = dd,
                     control = frmtmb_control(restarts = 8,
                                              grad_tol = 1e-12)))
fd <- capd$value
add(gc_row("cumulative n=6000, 8 restarts", 103, fd, capd$warnings))
cat("  objective delta after 8 restarts from the optimum:",
    format(fd$opt$objective - fc$opt$objective, digits = 6), "\n")
cat("  max|par delta|:",
    format(max(abs(fd$opt$par - fc$opt$par)), digits = 6), "\n")

## ------------------------------------------------------------------
hr("E. the same design at three sample sizes: does max|grad| track n?")
for (nn in c(200, 1000, 6000, 20000)) {
  set.seed(1030)
  d2 <- data.frame(x1 = rnorm(nn), x2 = rnorm(nn))
  e2 <- 0.8 * d2$x1 - 0.5 * d2$x2
  d2$yo <- cut(e2 + rlogis(nn), breaks = c(-Inf, -1, 0.5, 2, Inf),
               labels = FALSE)
  cp <- gc_catch(frm(bf(yo ~ x1 + x2), family = cumulative(), data = d2))
  add(gc_row(paste0("cumulative n = ", nn), 1030, cp$value, cp$warnings))
}

## ------------------------------------------------------------------
hr("F. a gaussian GLM: the closed form is exact")
for (nn in c(1000, 50000, 200000)) {
  set.seed(1040)
  d3 <- data.frame(x = rnorm(nn))
  d3$y <- rnorm(nn, 1 + 2 * d3$x, 2)
  cp <- gc_catch(frm(bf(y ~ x), family = gaussian(), data = d3))
  f3 <- cp$value
  add(gc_row(paste0("gaussian n = ", nn), 1040, f3, cp$warnings))
  lmf <- lm(y ~ x, data = d3)
  cat("  max relative coefficient gap against lm():",
      format(max(abs(coef(f3)[1:2] - coef(lmf)) / abs(coef(lmf))),
             digits = 4), "\n")
}

saveRDS(rows, file.path("dev", paste0("gradcheck-01-", which_lib, ".rds")))
cat("\nrows:", length(rows), "  warned:",
    sum(vapply(rows, function(r) isTRUE(r$warned), TRUE)), "\n")
