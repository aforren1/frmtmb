# What the user actually reads: the warning text on each branch, and the
# diagnose() block beside it.
#
#   Rscript dev/gradcheck-11-messages.R lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "lane"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n")

show <- function(lbl, expr) {
  cat("\n== ", lbl, "\n", sep = "")
  r <- gc_catch(expr)
  for (w in r$warnings) cat("  WARNING: ", w, "\n", sep = "")
  if (!length(r$warnings)) cat("  (no warnings)\n")
  if (!inherits(r$value, "gc_error")) {
    cat(paste0("  ", capture.output(diagnose(r$value))), sep = "\n")
  }
  invisible(NULL)
}

set.seed(101)
dd <- data.frame(x = rnorm(300))
dd$y <- rnorm(300, 1 + 2 * dd$x, 1)
show("an active bound (silent, but diagnose names it)",
     frm(bf(y ~ x), family = gaussian(), data = dd,
         prior = set_prior("", class = "b", ub = 0.1)))

set.seed(402)
n <- 1500
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rpois(n, exp(0.4 + 0.7 * d2$x - 0.4 * d2$z))
show("a loosened optimizer (warns, headroom measured)",
     frm(bf(y ~ x + z), family = poisson(), data = d2,
         control = frmtmb_control(
           restarts = 0,
           optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                          iter.max = 1000, eval.max = 1000))))

set.seed(201)
n <- 500
d3 <- data.frame(x = rnorm(n), z = rnorm(n))
d3$y <- rnorm(n, 3 + 2 * d3$x - 1.5 * d3$z, 1)
show("a far start (warns, curvature unusable)",
     frm(bf(y ~ x + z), family = gaussian(), data = d3,
         start = list(beta = c(50, -40, 30)),
         control = frmtmb_control(
           restarts = 0, optCtrl = list(iter.max = 1, eval.max = 3))))

set.seed(502)
n <- 1000
x1 <- rnorm(n)
d4 <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, 1e-5))
d4$y <- rnorm(n, 1 + 2 * d4$x1, 1)
show("a collinear ridge (warns on curvature, not size)",
     frm(bf(y ~ x1 + x2), family = gaussian(), data = d4,
         control = frmtmb_control(restarts = 0)))
