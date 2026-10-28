# Draws made without sample_prior on the lane's build against
# rellib-r14's, under fixed seeds: the md5 of each draws matrix, written
# to dev/priordraws-log/identity-<arm>-<rep>.txt. Run pinned to one core
# (dev/rtmb-pitfalls.md item 21).
#   Rscript dev/priordraws-identity.R <lane|base> <rep>
a <- commandArgs(TRUE)
arm <- a[1]
rep <- a[2]
.libPaths(c(if (identical(arm, "lane")) "C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
out <- file.path("C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/priordraws-log",
                 sprintf("identity-%s-%s.txt", arm, rep))
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(1)
d <- data.frame(x = rnorm(80), z = rnorm(80), g = factor(rep(1:8, each = 10)),
                m = sample(1:4, 80, TRUE), t = rep(1:10, 8))
d$y <- 1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(80)
runs <- list(
  re_fit = function() {
    fit <- q(frm(bf(y ~ x + (1 + x | g)), data = d))
    q(frm_sample(fit, chains = 2, iter = 400, seed = 3, refresh = 0))
  },
  mo_student = function() {
    q(frm_sample(bf(y ~ mo(m) + x), data = d, family = student(),
                 prior = set_prior("normal(0, 1)", class = "b"),
                 chains = 2, iter = 400, seed = 4, refresh = 0))
  },
  ar_cov = function() {
    q(frm_sample(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), data = d,
                 family = gaussian(), chains = 2, iter = 400, seed = 5,
                 refresh = 0))
  },
  laplace = function() {
    fit <- q(frm(bf(y ~ x + (1 | g)), data = d))
    q(frm_sample(fit, chains = 1, iter = 400, seed = 6, refresh = 0,
                 laplace = TRUE))
  })
lines <- c(sprintf("lib %s frmtmb.sample %s", find.package("frmtmb.sample"),
                   format(packageVersion("frmtmb.sample"))))
for (nm in names(runs)) {
  ds <- runs[[nm]]()
  m <- ds$draws
  f <- tempfile()
  saveRDS(unname(m), f, version = 2)
  lines <- c(lines, sprintf("%s %s %d x %d", nm, unname(tools::md5sum(f)),
                            nrow(m), ncol(m)))
}
writeLines(lines, out)
cat(lines, sep = "\n")
