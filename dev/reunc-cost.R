# Lane wt-reunc: what the group-effect draw costs predict().
#
# The arms are INTERLEAVED in one process and each block is grown past
# 1.2 s, because proc.time() ticks at 10 ms here (dev/lane-rules.md).
# The CONTROL is predict(re_formula = NA) on the same fit, which the
# lane does not change at all: it must report a ratio near 1 against the
# base build, and it says how much of any difference is the machine.
#
#   Rscript dev/reunc-cost.R <lib> <tag>
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]
tag <- a[2]
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
set.seed(11)
G <- 60
m <- 20
n <- G * m
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[d$g] + rnorm(n)
fit <- frm(bf(y ~ x + (1 | g)), data = d)
nd <- d[seq_len(50), c("x", "g")]
arms <- list(
  conditional = function() predict(fit, newdata = nd, ndraws = 200),
  population = function() predict(fit, newdata = nd, ndraws = 200,
                                  re_formula = NA),
  control = function() {   # fixed arithmetic, load-independent check
    s <- 0
    for (i in 1:120000) s <- s + sqrt(i)
    s
  })
best <- stats::setNames(rep(Inf, length(arms)), names(arms))
for (round in 1:5) {
  for (nm in names(arms)) {
    t0 <- proc.time()[["elapsed"]]
    reps <- 0L
    repeat {
      invisible(suppressWarnings(arms[[nm]]()))
      reps <- reps + 1L
      if (proc.time()[["elapsed"]] - t0 > 1.2) break
    }
    per <- (proc.time()[["elapsed"]] - t0) / reps
    best[[nm]] <- min(best[[nm]], per)
  }
}
cat(sprintf("%s (%s): %s\n", tag, find.package("frmtmb"),
            paste(sprintf("%s %.4f s", names(best), best), collapse = "; ")))
# the first predict() on a fit also pays for the joint precision, which
# an ML fit's sdreport does not carry; timed once, separately
f2 <- frm(bf(y ~ x + (1 | g)), data = d)
t0 <- proc.time()[["elapsed"]]
invisible(suppressWarnings(predict(f2, newdata = nd, ndraws = 10)))
cat(sprintf("%s: first predict() on a fresh fit %.3f s\n", tag,
            proc.time()[["elapsed"]] - t0))
