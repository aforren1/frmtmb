# Fit time of a by-split term against the same-size term without by.
#
# gr(g, by = f) adds nothing to the objective: it is K ordinary blocks
# where the term without by is one, so the Laplace Hessian over b keeps
# the same sparsity (block diagonal by level of g in both). What it adds
# is K - 1 covariance blocks' worth of theta. The arms are interleaved in
# one process, each timed over a block of repeated fits grown past 1.2 s,
# and the minimum of five rounds is kept; the control arm is the model
# without by timed a second time, so its ratio to the first must read 1.
#
# Run: Rscript dev/grby-timing.R > dev/grby-log/timing.txt
.libPaths(c("/opt/rlib/lane-grby", "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
set.seed(31)
ng <- 200
d <- data.frame(g = factor(rep(seq_len(ng), each = 10)))
d$f <- factor(ifelse(as.integer(d$g) <= ng / 2, "a", "b"))
d$x <- rnorm(nrow(d))
u <- cbind(rnorm(ng, 0, ifelse(seq_len(ng) <= ng / 2, 1, 0.4)),
           rnorm(ng, 0, 0.5))
d$y <- 1 + 0.5 * d$x + u[d$g, 1] + u[d$g, 2] * d$x + rnorm(nrow(d), 0, 0.5)
arms <- list(
  plain = function() frm(bf(y ~ x + (1 + x | g)) + gaussian(), data = d),
  by = function() frm(bf(y ~ x + (1 + x | gr(g, by = f))) + gaussian(),
                      data = d),
  control = function() frm(bf(y ~ x + (1 + x | g)) + gaussian(), data = d)
)
# one warm-up fit per arm first, so no arm pays the first taping, then
# one block size for every arm, grown on the slowest
for (a in arms) suppressWarnings(a())
k <- 1L
repeat {
  t <- min(vapply(arms, function(a) {
    system.time(for (i in seq_len(k)) suppressWarnings(a()))[["elapsed"]]
  }, 1))
  if (t > 1.2) break
  k <- k * 2L
}
reps <- vapply(arms, function(a) k, 1L)
per_fit <- matrix(NA_real_, 5, length(arms),
                  dimnames = list(NULL, names(arms)))
for (r in 1:5) {
  for (a in names(arms)) {
    per_fit[r, a] <- system.time(for (i in seq_len(reps[[a]])) {
      suppressWarnings(arms[[a]]())
    })[["elapsed"]] / reps[[a]]
  }
}
m <- apply(per_fit, 2, min)
cat("fits per timed block:", paste(names(reps), reps, collapse = ", "), "\n")
cat("seconds per fit, minimum of 5 rounds:\n")
print(round(m, 4))
cat("ratio by / plain", round(m[["by"]] / m[["plain"]], 3),
    "  control / plain", round(m[["control"]] / m[["plain"]], 3), "\n")
f_by <- arms$by()
f_pl <- arms$plain()
cat("theta length: by", length(f_by$estimates$theta), " plain",
    length(f_pl$estimates$theta), "\n")
# load-independent: how many times the optimizer asked for the objective
# and its gradient, each one inner Newton solve over b
cat("optimizer iterations: by", f_by$opt$iterations, " plain",
    f_pl$opt$iterations, "\n")
cat("objective / gradient evaluations: by",
    paste(f_by$opt$evaluations, collapse = " / "), " plain",
    paste(f_pl$opt$evaluations, collapse = " / "), "\n")
