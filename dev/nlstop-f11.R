# Lane nlstop, item 2584: F11 seed 8, bf(y ~ s(x1), sigma ~ s(x2)),
# gaussian. The default fit against starts with one or every log sd at
# -10 (the triage's check4.R), and where each better fit ends.
# Usage: NLSTOP_ARM=before|after Rscript dev/nlstop-f11.R [seeds]
source("dev/nlstop-common.R")
a <- commandArgs(TRUE)
seeds <- if (length(a)) eval(parse(text = a[1])) else 8
for (seed in seeds) {
  set.seed(seed)
  d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, verbose = FALSE))
  d$z <- runif(nrow(d)); d$g <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  d$g8 <- factor(sample(letters[1:8], nrow(d), TRUE))
  d$gr <- factor(rep(1:20, length.out = nrow(d)))
  d$y <- d$y + 0.6 * rnorm(20)[d$gr]
  fo <- bf(y ~ s(x1), sigma ~ s(x2))
  t0 <- proc.time()[["elapsed"]]
  r0 <- tt(frm(fo, data = d))
  t0 <- proc.time()[["elapsed"]] - t0
  th <- r0$estimates$theta
  cat(sprintf("seed %d default ll %.6f code %s t %.2f warn %d %s msg %d\n",
              seed, llv(r0), r0$opt$convergence, t0,
              length(attr(r0, "w")), wsum(r0), length(attr(r0, "m"))))
  cat("  theta:", format(th, digits = 5), "\n")
  for (w in attr(r0, "w")) cat("  W:", substr(gsub("\n", " ", w), 1, 160), "\n")
  best <- llv(r0)
  for (k in c(seq_along(th), 0L)) {
    st <- th
    if (k) st[k] <- -10 else st[] <- -10
    r <- tt(frm(fo, data = d, start = list(theta = st)))
    cat(sprintf("  start theta[%s] = -10: ll %.6f code %s warn %d theta %s\n",
                if (k) k else "all", llv(r), r$opt$convergence,
                length(attr(r, "w")),
                paste(format(r$estimates$theta, digits = 4), collapse = " ")))
    best <- max(best, llv(r), na.rm = TRUE)
  }
  cat(sprintf("seed %d gap %.6f\n", seed, llv(r0) - best))
}
