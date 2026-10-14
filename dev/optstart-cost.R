# Lane optstart: what the cold start costs on large designs, against the
# fit. make_start() on the assembled frame (minimum of 3 runs, elapsed
# seconds), predictor_units() alone (the solve the unit made before this
# lane), and one frm() fit. Run on the lane build and on rellib-r8.
#   Rscript dev/optstart-cost.R <lib or ""> <design> [fit]
# designs: wide (101000 by 50, the backlog case), long (1e6 by 5),
# factor (1e5 rows, a 500-level factor: a 1e5 by 501 design)
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
des <- a[2]
do_fit <- length(a) > 2 && a[3] == "fit"
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "design", des, "\n")
set.seed(1)
if (des == "wide") {
  n <- 101000
  p <- 49
  Z <- matrix(rnorm(n * p), n, p)
  colnames(Z) <- paste0("z", seq_len(p))
  d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
  d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
  fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                         "+ (1 | g)"))
} else if (des == "long") {
  n <- 1e6
  Z <- matrix(rnorm(n * 4), n, 4)
  colnames(Z) <- paste0("z", 1:4)
  d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
  d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
  fo <- y ~ z1 + z2 + z3 + z4 + (1 | g)
} else {
  n <- 1e5
  d <- data.frame(x = rnorm(n), f = factor(sample(500, n, TRUE)),
                  g = factor(rep(seq_len(100), length.out = n)))
  d$y <- 1 + 1e3 * d$x + rnorm(500)[d$f] + rnorm(100, 0, 0.5)[d$g] +
    rnorm(n)
  fo <- y ~ x + f + (1 | g)
}
bform <- ns$resolve_deferred_families(ns$expand_dot_bform(
  ns$as_bform(fo, gaussian()), d), d)
spec <- ns$parse_spec(bform)
frame <- ns$assemble_frame(spec, d, na.action = stats::na.omit,
                           data2 = list())
tm <- function(expr) {
  e <- substitute(expr)
  pf <- parent.frame()
  ts <- numeric(3)
  for (i in 1:3) ts[i] <- system.time(eval(e, pf))[["elapsed"]]
  min(ts)
}
t_ms <- tm(st <- ns$make_start(frame, NULL))
t_pu <- tm(u <- ns$predictor_units(frame))
cat(sprintf("make_start %.2f s | predictor_units %.2f s | unit %s\n",
            t_ms, t_pu, paste(signif(u, 3), collapse = ",")))
if (do_fit) {
  t_fit <- system.time(f <- suppressWarnings(frm(fo, data = d)))[["elapsed"]]
  cat(sprintf("fit %.1f s | code %d | objective %.10g\n", t_fit,
              f$opt$convergence, f$opt$objective))
}
