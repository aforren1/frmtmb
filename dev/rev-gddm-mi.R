# rev-gddm: the one place "design constant implies parameter constant"
# can be false. mi() replaces a missing covariate value with its own
# PARAMETER, one per missing entry, so two rows that are both NA in the
# model-frame column carry two different parameter values. The check
# treats two missing values as the same datum, by construction and on
# purpose, so it cannot see that.
#
# Seed 77. Arm from GDDM_LIB; default is this review's install.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, " frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n\n")

ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
set.seed(77)
n <- 60L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$cond <- rep(1:2, length.out = n)
# x is constant inside each condition where it is OBSERVED, and missing
# on four rows that are not their condition's first row
d$x <- ifelse(d$cond == 1L, 0.3, 0.7)
d$x[c(3L, 5L, 4L, 6L)] <- NA_real_
cat("x by condition, observed values:",
    paste(tapply(d$x, d$cond, function(z)
      paste(unique(z[!is.na(z)]), collapse = "/")), collapse = "  "),
    "\n")
cat("missing entries:", sum(is.na(d$x)), "on rows",
    paste(which(is.na(d$x)), collapse = ","), "\n\n")

sp <- bf(rt | vint(upper, cond) ~ mi(x), bs ~ 1, ndt ~ 1,
         bias = 0.5) + bf(x | mi() ~ 1, family = gaussian())
r <- tryCatch({
  frm(sp, family = gddm(control = ctl), data = d, dry_run = "frame")
  "ACCEPTED at frame assembly"
}, error = function(e) paste("refused:", conditionMessage(e)))
cat("mu ~ mi(x), four latent values inside two conditions:\n  ",
    substr(gsub("\\s+", " ", r), 1, 200), "\n", sep = "")

# the control: the same model with the missing entries observed and
# varying, which the check must refuse
d2 <- d; d2$x[c(3L, 5L)] <- 0.9; d2$x[c(4L, 6L)] <- 0.1
r2 <- tryCatch({
  frm(sp, family = gddm(control = ctl), data = d2, dry_run = "frame")
  "ACCEPTED"
}, error = function(e) paste("refused:",
                             substr(conditionMessage(e), 1, 70)))
cat("the same rows observed and varying:\n  ", r2, "\n", sep = "")

# the variant the NA-equality rule cannot see: EVERY row of a condition
# missing, so each row carries its own latent value and the column is
# all NA there
d3 <- d
d3$x[d3$cond == 1L] <- NA_real_
cat("\ncondition 1 entirely missing (", sum(is.na(d3$x)),
    " latent values in it):\n", sep = "")
r3 <- tryCatch({
  frm(sp, family = gddm(control = ctl), data = d3, dry_run = "frame")
  "ACCEPTED at frame assembly"
}, error = function(e) paste("refused:",
                            substr(conditionMessage(e), 1, 70)))
cat("  ", r3, "\n", sep = "")
