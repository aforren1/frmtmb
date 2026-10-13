# Reviewer of lane optrap, consolidation check: predictor_design_kappa()
# on ordinary designs with random effects, against the 1e5 that turns
# on newton_polish() for them. Frame only.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
kap <- frmtmb:::predictor_design_kappa
set.seed(1)
n <- 400
d <- data.frame(g = factor(rep(1:40, each = 10)),
                year = sample(1990:2020, n, TRUE), x = rnorm(n),
                f = factor(sample(letters[1:12], n, TRUE)),
                f2 = factor(sample(1:6, n, TRUE)),
                age = runif(n, 20, 80), xs = runif(n, 1, 2),
                big = 1e5 + runif(n), utime = 1.6e9 + runif(n, 0, 3e7))
d$y <- rnorm(n)
cases <- list(
  sleepstudy = list(Reaction ~ Days + (Days | Subject), sleepstudy),
  year = list(y ~ year + x + (1 | g), d),
  year_sq = list(y ~ year + I(year^2) + (1 | g), d),
  year_cub = list(y ~ year + I(year^2) + I(year^3) + (1 | g), d),
  age_sq = list(y ~ age + I(age^2) + (1 | g), d),
  age_cub = list(y ~ age + I(age^2) + I(age^3) + (1 | g), d),
  poly2_raw = list(y ~ xs + I(xs^2) + (1 | g), d),
  poly3_raw = list(y ~ xs + I(xs^2) + I(xs^3) + (1 | g), d),
  poly4_raw = list(y ~ xs + I(xs^2) + I(xs^3) + I(xs^4) + (1 | g), d),
  factors = list(y ~ f * f2 + (1 | g), d),
  big_offset = list(y ~ big + (1 | g), d),
  unix_time = list(y ~ utime + (1 | g), d),
  year_x_age = list(y ~ year * age + (1 | g), d))
for (nm in names(cases)) {
  fr <- frm(cases[[nm]][[1]], data = cases[[nm]][[2]], dry_run = "frame")
  k <- kap(fr)
  cat(sprintf("%-11s kappa %10.3g %s\n", nm, k,
              if (k > 1e5) "-> Newton steps run" else ""))
}
