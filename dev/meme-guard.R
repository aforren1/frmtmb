# Lane meme: models WITHOUT me() terms are unchanged. Run once per arm
# and diff the two outputs:
#   Rscript dev/meme-guard.R base > dev/meme-out/guard-base.txt
#   Rscript dev/meme-guard.R lane > dev/meme-out/guard-lane.txt
# Prints prior tables, fits at 17 digits, the objective and its
# gradient at a fixed point, and a digest of a seeded frm_sample() run.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
options(digits = 17)
set.seed(9)
n <- 80
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:8, 10)),
                t = rep(1:10, each = 8))
d$y <- 1 + 0.5 * d$x + rnorm(8, 0, 0.6)[d$g] + rnorm(n)
d$y2 <- d$y + 0.3 * d$z + rnorm(n)
d$sy <- 0.3
forms <- list(
  re = bf(y ~ x + (1 + x | g)),
  mv = bf(mvbind(y, y2) ~ x + (1 | g)) + set_rescor(TRUE),
  ar = bf(y ~ x + ar(time = t, gr = g, cov = TRUE)),
  mi = bf(y | mi(sy) ~ x))
pr <- set_prior("normal(0, 1)", class = "b") +
  set_prior("exponential(1)", class = "sd") +
  set_prior("lkj(2)", class = "cor")
for (nm in names(forms)) {
  f <- forms[[nm]]
  cat("==== ", nm, "\n")
  print(as.data.frame(default_prior(f, data = d)))
  print(as.data.frame(default_prior(f, data = d, route = "sample")))
  p <- if (identical(nm, "re")) pr else NULL
  fit <- suppressWarnings(frm(f, data = d, prior = p))
  print(unlist(fit$estimates))
  print(as.numeric(logLik(fit)))
  e <- fit$obj$env
  q <- e$last.par.best + 0.01
  print(e$f(q))
  print(as.numeric(e$f(q, order = 1)))
  if (!is.null(p)) print(as.data.frame(prior_summary(fit)))
}
ds <- suppressWarnings(suppressMessages(frm_sample(
  forms$re, data = d, prior = pr, chains = 1, iter = 300, seed = 3,
  refresh = 0, sample_prior = "yes")))
cat("draws digest", digest::digest(round(ds$draws, 10)), "\n")
print(colnames(ds$draws)[!grepl("^r_|^b\\[", colnames(ds$draws))])
print(as.data.frame(prior_summary(ds)))
print(names(prior_draws(ds)))
