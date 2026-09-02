# Probe: does the cluster-robust interval recover coverage where the
# model-based one loses it? Tunes the misspecification and the
# replicate count used by the seeded gate in test-sandwich.R.
#
# Run: Rscript dev/sandwich/probe-coverage.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
sink(file.path("dev", "sandwich", "probe-coverage.txt"), split = TRUE)

run <- function(label, nsim, G, m, slope_sd, x_between, het) {
  set.seed(17)
  n <- G * m
  gv <- factor(rep(seq_len(G), each = m))
  b <- 0.5
  cm <- cr <- logical(nsim)
  t0 <- proc.time()[["elapsed"]]
  for (s in seq_len(nsim)) {
    x <- rep(rnorm(G, 0, x_between), each = m) + rnorm(n, 0, 1)
    sl <- rep(rnorm(G, 0, slope_sd), each = m)
    sdg <- if (het) rep(runif(G, 0.2, 2.5), each = m) else 1
    y <- 1 + (b + sl) * x + rnorm(G, 0, 0.7)[gv] + rnorm(n, 0, sdg)
    dd <- data.frame(y = y, x = x, g = gv)
    ft <- try(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
                  REML = FALSE), silent = TRUE)
    if (inherits(ft, "try-error")) next
    ci_m <- confint(ft, parm = "x")
    V <- vcov_cluster(ft, ~ g, type = "CR1", full = TRUE)
    ci_r <- confint(ft, parm = "x", vcov = V)
    cm[s] <- ci_m[1, "lwr"] <= b && b <= ci_m[1, "upr"]
    cr[s] <- ci_r[1, "lwr"] <= b && b <= ci_r[1, "upr"]
  }
  cat(sprintf("%-34s model %.3f  robust %.3f  (%.1f s)\n",
              label, mean(cm), mean(cr),
              proc.time()[["elapsed"]] - t0))
}

run("het only", 60, 40, 6, 0.0, 0, TRUE)
run("random slope 0.6", 60, 40, 6, 0.6, 0, FALSE)
run("random slope 0.6 + het", 60, 40, 6, 0.6, 0, TRUE)
run("slope 0.8, x between 1.5, het", 60, 40, 6, 0.8, 1.5, TRUE)
run("slope 0.8, x between 1.5, het, m=12", 60, 40, 12, 0.8, 1.5, TRUE)
run("slope 1.0, x between 2, het, m=12", 60, 40, 12, 1.0, 2, TRUE)

sink()
