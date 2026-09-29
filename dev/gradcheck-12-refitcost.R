# The cost the new criterion adds to a REFIT LOOP.
#
# `frm_bootstrap()`, `influence()` and `frm_allfit()` wrap their refits in
# suppressWarnings(), which muffles the warning but does NOT stop it being
# built, so every replicate that trips the gradient trip-wire now pays for
# a Hessian. At n = 20000 that is 20 of 20 replicates on this design, so
# this is the worst case and not a typical one.
#
#   Rscript dev/gradcheck-12-refitcost.R base|lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "lane"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

run <- function(lbl, mk, nboot = 20L, rounds = 3L) {
  fit <- suppressWarnings(mk())
  d <- diagnose(fit, quiet = TRUE)
  best <- Inf
  for (k in seq_len(rounds)) {
    t0 <- proc.time()[["elapsed"]]
    invisible(suppressWarnings(frm_bootstrap(fit, nsim = nboot, seed = 1L)))
    best <- min(best, proc.time()[["elapsed"]] - t0)
  }
  # a control built from the same code that must report 1.0 against
  # itself: the plain refit loop, with no bootstrap machinery
  t1 <- proc.time()[["elapsed"]]
  for (k in seq_len(nboot)) invisible(suppressWarnings(mk()))
  ctl <- proc.time()[["elapsed"]] - t1
  cat(sprintf(paste0("%-30s np %2d  max|grad| %9.3e  bootstrap(%d)",
                     " %8.3f s  %d plain fits %8.3f s\n"),
              lbl, length(fit$opt$par), d$max_grad, nboot, best, nboot,
              ctl))
}

set.seed(1040)
n <- 20000
d1 <- data.frame(x = rnorm(n), z = rnorm(n))
d1$y <- rnorm(n, 1 + 2 * d1$x - 0.5 * d1$z, 2)
run("gaussian n=20000 (trips)",
    function() frm(bf(y ~ x + z), family = gaussian(), data = d1))

set.seed(1041)
n <- 400
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rnorm(n, 1 + 2 * d2$x - 0.5 * d2$z, 2)
run("gaussian n=400 (does not trip)",
    function() frm(bf(y ~ x + z), family = gaussian(), data = d2))

set.seed(1042)
n <- 20000
d3 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
e3 <- 0.8 * d3$x1 - 0.5 * d3$x2
d3$yo <- cut(e3 + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
             labels = FALSE)
run("cumulative n=20000 (trips)",
    function() frm(bf(yo ~ x1 + x2), family = cumulative(), data = d3),
    nboot = 10L, rounds = 2L)
