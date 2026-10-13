# brms 2.23.0's names for a logical and a character mo() multiplier.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(123)
lev <- c("a1", "a2", "a3", "a4")
income <- factor(sample(lev, 120, TRUE), levels = lev, ordered = TRUE)
dat <- data.frame(income = income, ls = c(30, 60, 70, 75)[income] + rnorm(120, sd = 7),
                  xc = sample(c("a", "b"), 120, TRUE))
dat$xl <- dat$xc == "a"
for (fo in list(ls ~ mo(income) * xl, ls ~ mo(income) * xc)) {
  dp <- default_prior(fo, data = dat)
  cat(deparse1(fo), ": b", dp$coef[dp$class == "b" & nzchar(dp$coef)],
      "| simo", dp$coef[dp$class == "simo"], "\n")
}
