# A yardstick for "the headroom predicts the realized drop", measured by
# the run instead of written down as 0.05.
#
# Along the Newton direction an exactly quadratic objective drops
# `D * (2t - t^2)`, so the drop at a half step is `0.75 * D`. The
# departure of the measured ratio from 0.75 is this objective's
# non-quadraticity AT THIS POINT, and the headroom's own error is the
# same third-order term. So the assertion can be a ratio between two
# things the run measures.
#
#   Rscript dev/gradcheck-14-quadyard.R lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "lane"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

probe <- function(lbl, fit) {
  p <- fit$opt$par
  g <- drop(fit$obj$gr(p))
  H <- fit$obj$he(p)
  ch <- chol((H + t(H)) / 2)
  step <- -backsolve(ch, backsolve(ch, g, transpose = TRUE))
  f0 <- as.numeric(fit$obj$fn(p))
  d1 <- f0 - as.numeric(fit$obj$fn(p + step))
  d05 <- f0 - as.numeric(fit$obj$fn(p + 0.5 * step))
  head <- suppressWarnings(diagnose(fit, quiet = TRUE)$grad_headroom)
  # non-quadraticity, measured: 0.75 is the exact-quadratic value
  nonq <- abs(d05 / d1 / 0.75 - 1)
  err <- abs(head / d1 - 1)
  cat(sprintf(paste0("%-34s d1 %11.5g  d05/d1 %8.5f  nonquad %9.3e",
                     "  err %9.3e  err/nonquad %8.3f\n"),
              lbl, d1, d05 / d1, nonq, err, err / nonq))
  invisible(NULL)
}

loose <- function(rt) {
  frmtmb_control(restarts = 0,
                 optCtrl = list(rel.tol = rt, x.tol = rt,
                                iter.max = 2000, eval.max = 2000))
}

set.seed(402)
n <- 1500
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rpois(n, exp(0.4 + 0.7 * d2$x - 0.4 * d2$z))
for (rt in c(1e-2, 1e-3, 1e-4)) {
  probe(paste0("poisson rel.tol = ", rt),
        suppressWarnings(frm(bf(y ~ x + z), family = poisson(), data = d2,
                             control = loose(rt))))
}

set.seed(9301)
n <- 3000
d3 <- data.frame(x = rnorm(n), z = rnorm(n))
d3$y <- rpois(n, exp(0.4 + 0.9 * d3$x - 0.6 * d3$z))
probe("poisson bounded, rel.tol = 1e-2",
      suppressWarnings(frm(bf(y ~ x + z), family = poisson(), data = d3,
                           prior = set_prior("", class = "b", coef = "x",
                                             ub = 0.1),
                           control = loose(1e-2))))

set.seed(204)
n <- 6000
d4 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
e4 <- 0.8 * d4$x1 - 0.5 * d4$x2
d4$yo <- cut(e4 + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
             labels = FALSE)
for (cap in c(4L, 8L)) {
  probe(paste0("cumulative n=6000, cap ", cap),
        suppressWarnings(frm(
          bf(yo ~ x1 + x2), family = cumulative(), data = d4,
          control = frmtmb_control(
            restarts = 0,
            optCtrl = list(iter.max = cap, eval.max = cap * 3)))))
}
