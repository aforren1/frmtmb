# One row of gradcheck-08 disagreed: at collinearity 1e-6 the headroom
# predicts 0.528 and a harder optimization from the same point realized
# 2.3e-8. Either the Hessian there is too ill-conditioned to invert, or
# nlminb cannot walk the direction. Settled by taking the step and
# evaluating the objective, which needs no optimizer at all.
#
#   Rscript dev/gradcheck-09-newtonstep.R lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "lane"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

step_probe <- function(lbl, fit) {
  p <- fit$opt$par
  g <- drop(fit$obj$gr(p))
  H <- tryCatch(fit$obj$he(p), error = function(e) NULL)
  src <- "he()"
  if (is.null(H)) {
    H <- stats::optimHess(p, function(q) as.numeric(fit$obj$fn(q)),
                          function(q) drop(fit$obj$gr(q)))
    src <- "optimHess"
  }
  H <- (H + t(H)) / 2
  ch <- chol(H)
  d <- -backsolve(ch, backsolve(ch, g, transpose = TRUE))
  pred <- 0.5 * sum(backsolve(ch, g, transpose = TRUE)^2)
  f0 <- as.numeric(fit$obj$fn(p))
  best <- f0
  bt <- 0
  for (t in c(1, 0.5, 0.25, 0.1, 2, 4)) {
    ft <- tryCatch(as.numeric(fit$obj$fn(p + t * d)),
                   error = function(e) NA_real_)
    if (is.finite(ft) && ft < best) {
      best <- ft
      bt <- t
    }
  }
  cat(sprintf("%-26s H from %-10s rcond %9.3e  |step| %10.3e\n",
              lbl, src, 1 / kappa(H), sqrt(sum(d^2))))
  cat(sprintf(paste0("%-26s predicted %12.6g  best step drop %12.6g",
                     " at t = %g  ratio %9.4f\n"),
              "", pred, f0 - best, bt,
              (f0 - best) / pred))
  invisible(NULL)
}

for (eps in c(1e-5, 1e-6, 1e-7)) {
  set.seed(502)
  n <- 1000
  x1 <- rnorm(n)
  d <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, eps))
  d$y <- rnorm(n, 1 + 2 * d$x1, 1)
  f <- suppressWarnings(frm(bf(y ~ x1 + x2), family = gaussian(), data = d,
                            control = frmtmb_control(restarts = 0)))
  step_probe(paste0("collinear eps = ", eps), f)
}

# the well-conditioned controls, where gradcheck-08 already agreed
set.seed(402)
n <- 1500
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rpois(n, exp(0.4 + 0.7 * d2$x - 0.4 * d2$z))
for (rt in c(1e-2, 1e-3)) {
  f <- suppressWarnings(frm(
    bf(y ~ x + z), family = poisson(), data = d2,
    control = frmtmb_control(restarts = 0,
                             optCtrl = list(rel.tol = rt, x.tol = rt,
                                            iter.max = 1000,
                                            eval.max = 1000))))
  step_probe(paste0("poisson rel.tol = ", rt), f)
}
