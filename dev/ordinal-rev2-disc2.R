# Reviewer re-check: allFit output without its timing column (base vs
# lane); emmeans on a modeled disc; coef() of a user-fixed dpar on base
# (student nu = 4) for the convention. Data seed 20261013.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
b <- readRDS(file.path(wt, "dev/ordinal-rev2-disc-base.rds"))$allfit
l <- readRDS(file.path(wt, "dev/ordinal-rev2-disc-lane.rds"))$allfit
strip <- function(x) sub("[0-9.]+$", "", trimws(x))
cat("allFit print without the seconds column identical:",
    identical(strip(b), strip(l)), "\n")
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib", "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(emmeans)})
set.seed(20261013)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(letters[1:6], n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x +
  rnorm(6, 0, 0.4)[as.integer(d$g)]
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
f <- frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(), data = d)
print(emmeans(f, ~ z, dpar = "disc", at = list(z = c(-1, 1))))
d$yg <- rnorm(n) + d$x
fs <- frm(bf(yg ~ x, nu = 4), family = student(), data = d)
cat("student nu = 4: coef names:", names(coef(fs)), "| flatten:",
    names(fixef(fs, flatten = TRUE)), "| variables:", variables(fs), "\n")
