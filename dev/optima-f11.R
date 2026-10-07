# Lane optima, item 4b: fixes' F11 seed 8, bf(y ~ s(x1), sigma ~
# s(x2)) on gamSim eg 6 (the builder of dev/fixes-p2-f11.R). The
# default fit, the fit from the variance-component boundary (each
# smoothing SD's log started at -10, alone and together), and a
# multi-start of 20 starts jittered about the default's estimates
# (theta by N(0, 2)), against which the default is judged.
#   Rscript dev/optima-f11.R base|lane [seed]
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seed <- if (length(args) > 1) as.integer(args[2]) else 8L
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
gs <- function(eg, seed, n = 200) {
  set.seed(seed)
  d <- suppressMessages(mgcv::gamSim(eg = eg, n = n, verbose = FALSE))
  d$z <- runif(nrow(d))
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  d$g8 <- factor(sample(letters[1:8], nrow(d), TRUE))
  d$gr <- factor(rep(1:20, length.out = nrow(d)))
  d$y <- d$y + 0.6 * stats::rnorm(20)[d$gr]
  d
}
d <- gs(6, seed)
fo <- bf(y ~ s(x1), sigma ~ s(x2))
one <- function(lab, start = NULL) {
  f <- suppressWarnings(frm(fo, data = d, start = start))
  cat(sprintf("%-28s conv %d logLik %.6f theta %s evals %s\n", lab,
              f$opt$convergence, as.numeric(logLik(f)),
              paste(sprintf("%.2f", f$estimates$theta), collapse = " "),
              f$opt$evals))
  f
}
f0 <- one("default")
th <- f0$estimates$theta
for (k in seq_along(th)) {
  st <- th
  st[k] <- -10
  one(sprintf("theta[%d] from -10", k), list(theta = st))
}
one("both theta from -10", list(theta = rep(-10, length(th))))
set.seed(1)
ll <- vapply(1:20, function(i) {
  st <- list(theta = th + stats::rnorm(length(th), 0, 2))
  f <- suppressWarnings(frm(fo, data = d, start = st))
  as.numeric(logLik(f))
}, 0)
cat("multi-start (20): max", format(max(ll), digits = 10), "; default",
    format(as.numeric(logLik(f0)), digits = 10), "; gap",
    format(max(ll) - as.numeric(logLik(f0)), digits = 4), "\n")
print(table(round(ll, 4)))
