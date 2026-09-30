# Reviewer, lane ordinal: show that the mutants the suite misses are
# behaviorally wrong. Usage: Rscript ... <lane|M2_acat_logit_nodisc|M3b_stz_inverse>
# Data seed 62 (the sample test's), sampler seed 3; core data 20260930.
# Output: dev/ordinal-rev-log-mutdemo-<arm>.txt
arm <- commandArgs(TRUE)[1]
sp <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad"
libs <- c("C:/Users/adf44/source/r/wt-ordinal-lib",
          "C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm != "lane") libs <- c(file.path(sp, "ordrev-mutlib", arm), libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
# the core test's data and model (test-ordinal-disc-thres.R, block 1)
set.seed(20260930)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
f <- frm(bf(y ~ x, disc ~ 0 + z), family = acat(), data = d)
cat("acat logit disc ~ 0 + z: logLik", format(as.numeric(logLik(f)), digits = 12),
    " disc coefficient", fixef(f, flatten = TRUE)[["disc_z"]], "\n")
f2 <- frm(bf(y ~ x, disc ~ 0 + z), family = cratio(), data = d)
cat("cratio logit disc ~ 0 + z: logLik", format(as.numeric(logLik(f2)), digits = 12),
    " disc coefficient", fixef(f2, flatten = TRUE)[["disc_z"]], "\n")
# the sample test's sum-to-zero draws
set.seed(62)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
fit <- frm(y ~ x, family = sratio(threshold = "sum_to_zero"), data = d)
ds <- suppressWarnings(suppressMessages(frm_sample(fit, chains = 1, iter = 300,
                                                   refresh = 0, seed = 3)))
ep <- posterior_epred(ds, ndraws = 3)
m <- as.matrix(ds, variable = c(paste0("b_Intercept[", 1:4, "]"), "b_x"))
# brms's own density at draw 1's stored thresholds
P <- brms:::dsratio(1:5, eta = d$x * m[1, "b_x"],
                    thres = matrix(m[1, 1:4], n, 4, byrow = TRUE))
cat("sratio stz draws: max |epred(draw 1) - brms at stored columns|",
    max(abs(ep[1, , ] - P)), "\n")
ll <- log_lik(ds)
cat("log_lik[1, 1:3]:", ll[1, 1:3], "\n")
