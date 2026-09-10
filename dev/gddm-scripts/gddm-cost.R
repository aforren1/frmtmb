# lane gddm: what the alternative to refusing would cost.
#
# The alternative to refusing a dpar that varies within a condition is
# to make the density READ it, and there is only one way to do that: one
# Fokker-Planck solve per distinct parameter vector. The tape cannot
# compare parameter values, so the distinct vectors have to be found
# from the DESIGN at frame assembly, which is the operation
# gddm_conditions() already performs and which the user can perform
# themselves. So "read it" IS "refine the condition index", and its cost
# is exactly one solve per condition.
#
# The exact quantity is load-independent and needs no clock: solves per
# likelihood evaluation IS the number of conditions. The clock is here
# to say what one solve costs.
#
# Two instrument problems this script hit, both recorded because both
# produced a confident wrong number first:
#
#  1. obj$fn(p) at an UNCHANGED p returns TMB's cached value. The first
#     version measured 10 us for 120 Fokker-Planck solves. Every call
#     now moves one coordinate.
#  2. the arms were run in one order and the control arm always followed
#     the heaviest one, so it carried its garbage collection: the
#     control came back at 1.9 where 1.0 is right. The order is now a
#     palindrome and gc() runs before each timed arm, so the two ends
#     see the same conditions.
#
# Seed 77. Arm chosen by GDDM_LIB.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "  frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

set.seed(77)
n <- 120L
grid <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))

mk <- function(ncond) {
  dd <- d
  dd$cond <- ((seq_len(n) - 1L) %% ncond) + 1L
  dd$k <- factor(dd$cond)
  frm(bf(rt | vint(upper, cond) ~ 0 + k, bs ~ 1, ndt ~ 1, bias = 0.5),
      family = gddm(control = grid), data = dd,
      dry_run = "objective")$obj
}

# A block grows until it takes at least 1.2 s, so the 10 ms clock tick
# on this machine cannot be the measurement.
block <- function(f, target = 1.2) {
  reps <- 1L
  repeat {
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(reps)) f(i)
    el <- proc.time()[["elapsed"]] - t0
    if (el >= target) return(el / reps)
    reps <- reps * 2L
  }
}

conds <- c(2L, 4L, 8L, 20L, 60L, 120L)
arms <- c(conds, rev(conds)[-1L])     # palindrome: 2 .. 120 .. 2
rounds <- 3L
bt <- rep(Inf, length(arms))
et <- rep(Inf, length(arms))
for (r in seq_len(rounds)) {
  for (a in seq_along(arms)) {
    invisible(gc(FALSE))
    t0 <- proc.time()[["elapsed"]]
    ob <- mk(arms[[a]])
    b <- proc.time()[["elapsed"]] - t0
    p <- ob$par
    j <- which(names(p) == "betad")[[1L]]
    # a moved coordinate, because TMB returns the cached value at an
    # unchanged parameter vector
    e <- block(function(i) {
      q <- p
      q[[j]] <- p[[j]] + i * 1e-9
      ob$fn(q)
    })
    if (b < bt[[a]]) bt[[a]] <- b
    if (e < et[[a]]) et[[a]] <- e
  }
}

# The table prints the ASCENDING arm; the descending half is read only
# through the controls below, which compare the two ends. An earlier
# version computed a fold it never used, which the review caught.

cat("\n== 120 rows, grid dt = 0.05, ny = 51, t_max = 2\n")
cat(sprintf("%8s %13s %13s %13s %13s\n", "ncond", "tape build s",
            "one fn s", "ms per solve", "fn vs ncond=2"))
for (a in seq_along(conds)) {
  cat(sprintf("%8d %13.3f %13.4f %13.3f %13.1f\n", conds[[a]],
              bt[[a]], et[[a]], 1e3 * et[[a]] / conds[[a]],
              et[[a]] / et[[1L]]))
}
cat(sprintf("\ncontrol, ncond = 2 at both ends of the palindrome: fn %.3f, tape %.3f (1.0 is right)\n",
            et[[length(arms)]] / et[[1L]], bt[[length(arms)]] / bt[[1L]]))
cat(sprintf("control, ncond = 4 at both ends: fn %.3f, tape %.3f\n",
            et[[length(arms) - 1L]] / et[[2L]],
            bt[[length(arms) - 1L]] / bt[[2L]]))
cat(sprintf("\none solve per ROW against one per condition, 120 rows: fn %.1fx, tape build %.1fx\n",
            et[[length(conds)]] / et[[1L]],
            bt[[length(conds)]] / bt[[1L]]))
cat("solves per likelihood evaluation IS ncond, exactly, no clock\n")
