# Reviewer of lane illse, round 2: seed 1 and 4 of dev/illse-rev2-flat.R,
# start log sd -15, iter.max 3: the warnings in full and the gradient.
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
for (seed in c(1, 4)) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.4)
  r <- cap(frm(y ~ x + (1 | g), data = d, start = list(theta = -15),
               control = frmtmb_control(optCtrl = list(iter.max = 3, eval.max = 1000))))
  f <- r$v
  cat("seed", seed, "code", f$opt$convergence, f$opt$message, "\n")
  print(f$opt$par)
  print(setNames(as.numeric(f$obj$gr(f$opt$par)), names(f$opt$par)))
  cat(paste0(" W: ", r$w), sep = "\n")
  print(sqrt(diag(suppressWarnings(vcov(f, full = TRUE)))))
}
