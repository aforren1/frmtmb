# Reviewer of lane seflags, final check: the fit warnings of
# dev/seflags-rev3-illid.R's no-random-effect fits (seeds 4, 5), where
# the op3 build gives NaN SEs everywhere with no warning on predictions.
#   Rscript dev/seflags-rev3-op3w.R <arm>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
for (s in 4:5) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                  z = stats::rnorm(200),
                  w = factor(sample(c("a", "b"), 200, TRUE)))
  d$y <- sin(3 * d$x) + 0.3 * d$z + 0.5 * (d$w == "b") +
    stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  r <- cap(frm(stats::as.formula(paste("y ~", paste0("x", 1:8, collapse = " + "),
                                       "+ z + w")), data = d))
  cat("seed", s, "code", r$v$opt$convergence, "lost", lost_of(r$v), "\n")
  for (w in r$w) cat("  W:", substr(w, 1, 300), "\n")
}
