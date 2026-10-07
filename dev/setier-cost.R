# Lane setier: what the tier-3-first rule costs a fit. One arm per
# process (base and lane cannot share one); dev/setier-cost.sh
# alternates the arms. Each model is fitted in a block of repeats grown
# past 1.2 s, and the per-fit time of the block is printed.
#   Rscript dev/setier-cost.R <lib or "base"> <round>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
arm <- if (identical(args[1], "base")) "base" else "lane"
d <- lme4::sleepstudy
d$a <- factor(d$Days %% 3)
set.seed(2)
dg <- data.frame(g = factor(rep(1:30, each = 8)), f = factor(rep(1:30, 8)),
                 x = rnorm(240))
dg$y <- rbinom(240, 1, plogis(0.5 * dg$x + rnorm(30)[dg$g]))
set.seed(3)
dn <- data.frame(x = rnorm(300), z = rnorm(300))
dn$y <- rpois(300, exp(0.3 + 0.5 * dn$x - 0.2 * dn$z))
models <- list(
  # a boundary fit: the prescreen, tier 3 and the at-zero probe run
  nested = function() frm(Reaction ~ Days + (1 | Subject/a), data = d),
  # healthy, fit-time check: the prescreen's eigenvalues only
  slope = function() frm(Reaction ~ Days + (Days | Subject), data = d),
  # 32 outer parameters, the check waits; summary() builds the Hessian
  glmm_f30 = function() {
    f <- frm(y ~ x + f + (1 | g), family = bernoulli(), data = dg)
    summary(f)
    f
  },
  # control: no random effects, nothing in the lane's path changes
  control = function() frm(y ~ x + z, family = poisson(), data = dn))
for (nm in names(models)) {
  f <- models[[nm]]
  suppressMessages(suppressWarnings(f()))
  reps <- 1L
  repeat {
    t <- system.time(for (i in seq_len(reps)) {
      suppressMessages(suppressWarnings(f()))
    })[["elapsed"]]
    if (t >= 1.2) break
    reps <- reps * 2L
  }
  cat(sprintf("COST %s %s round %s per_fit %.5f reps %d\n", arm, nm,
              args[2], t / reps, reps))
}
