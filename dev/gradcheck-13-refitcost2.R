# The refit-loop cost, remeasured the way the instrument rules require.
#
# dev/gradcheck-12-refitcost.R timed two BUILDS in two processes, which
# confounds the reading with everything else that differs between them,
# and offered an n = 400 arm as a control that is 13 ticks of a 10 ms
# clock. This times the two READINGS inside ONE process on ONE build:
# `restarts = 0` in both arms so `grad_tol` governs only the warning
# path, the "off" arm sets `grad_tol` to ten times the fit's own gradient
# so stage 2 never runs, the arms are interleaved per round, blocks are
# grown past 1.2 s, the minimum of 5 rounds is taken, and a SECOND "off"
# arm is the control that must report 1.0 on identical work.
#
#   Rscript dev/gradcheck-13-refitcost2.R lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "lane"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

run <- function(lbl, mk, nboot = 20L, rounds = 5L, min_block = 1.2) {
  fit <- suppressWarnings(mk(1e-3))
  g <- suppressWarnings(diagnose(fit, quiet = TRUE)$max_grad)
  # `on` is the default reading; `off` puts grad_tol above the fit's own
  # gradient, so the trip-wire never fires and no Hessian is built
  arms <- list(on = 1e-3, off = 10 * g, ctl = 10 * g)
  best <- c(on = Inf, off = Inf, ctl = Inf)
  for (k in seq_len(rounds)) {
    for (a in names(arms)) {
      f <- suppressWarnings(mk(arms[[a]]))
      t0 <- proc.time()[["elapsed"]]
      reps <- 0L
      repeat {
        invisible(suppressWarnings(frm_bootstrap(f, nsim = nboot,
                                                 seed = 1L)))
        reps <- reps + 1L
        if (proc.time()[["elapsed"]] - t0 > min_block) break
      }
      best[[a]] <- min(best[[a]],
                       (proc.time()[["elapsed"]] - t0) / reps)
    }
  }
  cat(sprintf(paste0("%-30s max|grad| %9.3e  on %8.3f s  off %8.3f s",
                     "  ratio %6.3f  control %6.3f\n"),
              lbl, g, best[["on"]], best[["off"]],
              best[["on"]] / best[["off"]],
              best[["ctl"]] / best[["off"]]))
  invisible(NULL)
}

set.seed(1040)
n <- 20000
d1 <- data.frame(x = rnorm(n), z = rnorm(n))
d1$y <- rnorm(n, 1 + 2 * d1$x - 0.5 * d1$z, 2)
run("gaussian n=20000 (trips)", function(gt) {
  frm(bf(y ~ x + z), family = gaussian(), data = d1,
      control = frmtmb_control(restarts = 0, grad_tol = gt))
})

set.seed(1042)
n <- 20000
d3 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
e3 <- 0.8 * d3$x1 - 0.5 * d3$x2
d3$yo <- cut(e3 + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
             labels = FALSE)
run("cumulative n=20000 (trips)", function(gt) {
  frm(bf(yo ~ x1 + x2), family = cumulative(), data = d3,
      control = frmtmb_control(restarts = 0, grad_tol = gt))
}, nboot = 10L, rounds = 3L)

# a design that does NOT trip: the two readings run the same code, so the
# effect arm must report 1.0 as well. Kept large enough to clear the
# clock, which the n = 400 arm of gradcheck-12 did not.
set.seed(1043)
n <- 6000
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rnorm(n, 1 + 2 * d2$x - 0.5 * d2$z, 2)
run("gaussian n=6000, grad_tol raised", function(gt) {
  frm(bf(y ~ x + z), family = gaussian(), data = d2,
      control = frmtmb_control(restarts = 0, grad_tol = max(gt, 1)))
})
