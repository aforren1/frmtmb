# Is the headroom the criterion warns on REAL? The predicted gain from
# one Newton step, against the gain a harder optimization actually gets
# from the same stopping point.
#
# The ill-conditioned rows matter most: there the Hessian is nearly
# singular, so `g' H^-1 g` could be an artifact of inverting it.
#
#   Rscript dev/gradcheck-08-calib.R lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "lane"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

cat(sprintf("%-32s %13s %13s %9s\n", "case", "predicted", "realized",
            "ratio"))
row <- function(lbl, fit, harder) {
  pr <- gc_probe(fit)
  real <- fit$opt$objective - harder$opt$objective
  cat(sprintf("%-32s %13.6g %13.6g %9.3f\n", lbl, pr$decr, real,
              pr$decr / real))
  invisible(NULL)
}

hard <- frmtmb_control(restarts = 12, grad_tol = 1e-14,
                       optCtrl = list(iter.max = 20000,
                                      eval.max = 20000,
                                      rel.tol = 1e-14, x.tol = 1e-14))

## the ill-conditioned ridge, where the decrement could be an artifact
for (eps in c(1e-5, 1e-6)) {
  set.seed(502)
  n <- 1000
  x1 <- rnorm(n)
  d <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, eps))
  d$y <- rnorm(n, 1 + 2 * d$x1, 1)
  f <- suppressWarnings(frm(bf(y ~ x1 + x2), family = gaussian(), data = d,
                            control = frmtmb_control(restarts = 0)))
  hf <- suppressWarnings(frm(bf(y ~ x1 + x2), family = gaussian(),
                            data = d, control = hard,
                            start = list(beta = unname(f$opt$par[1:3]))))
  row(paste0("collinear eps = ", eps), f, hf)
}

## a loosened nlminb on a GLM and on a Laplace GLMM
set.seed(402)
n <- 1500
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rpois(n, exp(0.4 + 0.7 * d2$x - 0.4 * d2$z))
for (rt in c(1e-2, 1e-3, 1e-5)) {
  f <- suppressWarnings(frm(
    bf(y ~ x + z), family = poisson(), data = d2,
    control = frmtmb_control(restarts = 0,
                             optCtrl = list(rel.tol = rt, x.tol = rt,
                                            iter.max = 1000,
                                            eval.max = 1000))))
  hf <- suppressWarnings(frm(bf(y ~ x + z), family = poisson(), data = d2,
                             control = hard,
                             start = list(beta = unname(f$opt$par[1:3]))))
  row(paste0("poisson rel.tol = ", rt), f, hf)
}

set.seed(403)
n <- 1200
d3 <- data.frame(x = rnorm(n), g = factor(rep(1:40, 30)))
d3$y <- rpois(n, exp(0.4 + 0.6 * d3$x + rnorm(40, 0, 0.7)[d3$g]))
for (rt in c(1e-2, 1e-3, 1e-4)) {
  f <- suppressWarnings(frm(
    bf(y ~ x + (1 | g)), family = poisson(), data = d3,
    control = frmtmb_control(restarts = 0,
                             optCtrl = list(rel.tol = rt, x.tol = rt,
                                            iter.max = 1000,
                                            eval.max = 1000))))
  hf <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = poisson(),
                             data = d3, control = hard,
                             start = list(beta = unname(f$opt$par[1:2]),
                                          theta = unname(f$opt$par[3]))))
  row(paste0("poisson GLMM rel.tol = ", rt), f, hf)
}

## a CORRECT fit that trips the trip-wire: the realized gain has to be
## as small as the predicted one, or the criterion is letting a real
## shortfall through
set.seed(1040)
for (nn in c(50000L, 200000L)) {
  d <- data.frame(x = rnorm(nn))
  d$y <- rnorm(nn, 1 + 2 * d$x, 2)
  f <- suppressWarnings(frm(bf(y ~ x), family = gaussian(), data = d))
  hf <- suppressWarnings(frm(bf(y ~ x), family = gaussian(), data = d,
                             control = hard,
                             start = list(beta = unname(f$opt$par[1:2]))))
  row(paste0("gaussian n = ", nn, " (correct)"), f, hf)
}
