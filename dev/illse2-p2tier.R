# Lane illse2, item 2: test-se-tier.R's "an sd the fit left near zero
# short of the maximum" fits, from a group sd of exp(-20), against lme4.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-p2tier.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
for (seed in c(36, 19)) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  r <- cap(frm(y ~ x + (1 | g), data = d, start = list(theta = -20)))
  f <- r$v
  l <- lme4::lmer(y ~ x + (1 | g), data = d, REML = FALSE)
  cat(sprintf(paste0("seed %d | code %d | logLik - lme4 %+.2e | sd %.4g ",
                     "(lme4 %.4g) | escape %s | lost %s | warn %s\n"),
              seed, f$opt$convergence,
              as.numeric(logLik(f)) - as.numeric(logLik(l)),
              exp(f$opt$par[names(f$opt$par) == "theta"]),
              attr(lme4::VarCorr(l)$g, "stddev"),
              if (is.null(f$opt$sd_escape)) "-" else
                sprintf("%.3g", f$opt$sd_escape[["gain"]]),
              lost_of(f), wtag(r$w)))
}
