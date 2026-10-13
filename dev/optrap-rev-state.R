# Reviewer of lane optrap, item 4: after a fit, are the objective's last
# and best points the reported optimum, and are the random-effect modes
# the ones at the optimum? `ls ~ mo(income) * age + (1 | grp)`, the
# mo-study design with 10 groups added, seeds 1 to 60, per build.
#   Rscript dev/optrap-rev-state.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
n_last <- n_best <- n_modes <- n_search <- 0L
for (s in 1:60) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  d <- data.frame(income, grp = factor(rep(1:10, 10)))
  d$ls <- c(30, 60, 70, 75)[income] + rnorm(10, 0, 3)[d$grp] +
    rnorm(100, sd = 7)
  d$age <- rnorm(100, mean = 40, sd = 10)
  f <- suppressWarnings(frm(ls ~ mo(income) * age + (1 | grp), data = d))
  e <- f$obj$env
  lf <- e$lfixed()
  p <- unname(f$opt$par)
  at_last <- identical(unname(e$last.par[lf]), p)
  at_best <- identical(unname(e$last.par.best[lf]), p)
  r_now <- e$last.par.best[-lf]
  f$obj$fn(f$opt$par)
  r_opt <- f$obj$env$last.par[-lf]
  modes_ok <- max(abs(r_now - r_opt)) < 1e-8
  n_last <- n_last + at_last
  n_best <- n_best + at_best
  n_modes <- n_modes + modes_ok
  n_search <- n_search + !is.null(f$opt$mo_search)
  if (!(at_last && at_best && modes_ok)) {
    cat(sprintf("seed %d: last %s best %s modes %s (max diff %.3g)\n", s,
                at_last, at_best, modes_ok, max(abs(r_now - r_opt))))
  }
}
cat(a[2], "of 60: last point at the optimum", n_last, "| best point",
    n_best, "| modes those at the optimum", n_modes,
    "| fits whose mo search kept a run", n_search, "\n")
