# What more draws measure on the RL model's 100-trial shifts, and the worst
# group's ESS of the 30 x 2 design over data seeds. Log dev/impsamp-out/rl2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
sys.source(system.file("rl", "rw-delta.R", package = "frmtmb",
                       mustWork = TRUE), envir = globalenv())
sq <- function(e) suppressWarnings(suppressMessages(e))
form <- bf(choice | reward(pay1, pay2) ~ condition + (1 | p | id),
           beta ~ 1 + (1 | p | id))
for (s in c(4201L, 4205L)) {
  set.seed(s)
  dd <- rl_simulate(rl_bandit_design(40, 100), seed = s)[[1L]]
  fam <- rw_delta(subject = id, trial = trial)
  lap <- sq(frm(form, family = fam, data = dd))
  cat(sprintf("seed %d Laplace log sd %.4f\n", s,
              log(VarCorr(lap)$id$sd[1, "Estimate"])))
  for (N in c(100L, 400L, 1600L)) {
    f <- sq(frm(form, family = fam, data = dd, importance = N))
    cat(sprintf("  N=%d log sd %.4f mcse %.4f essmin %.3f rounds %d\n", N,
                log(VarCorr(f)$id$sd[1, "Estimate"]), f$importance$mcse,
                f$importance$ess_min, f$importance$rounds))
  }
}
ess <- numeric(0)
for (s in 1:10) {
  set.seed(s + 4L)
  ng <- 30L
  g <- factor(rep(seq_len(ng), each = 2L))
  x <- rnorm(2L * ng)
  du <- data.frame(y = rbinom(2L * ng, 1, 0.4), x = x, g = g)
  f <- sq(frm(bf(y ~ x + (x | g)) + bernoulli(), data = du, importance = 500L))
  ess <- c(ess, f$importance$ess_min)
}
cat("30 x 2, data seeds 5:14, 500 draws, worst group:",
    paste(format(ess, digits = 2), collapse = " "), "\n")
