# Reviewer of lane sampleprior: the lane's posterior-agreement cases
# (dev/sampleprior-agree-common.R, read only) rerun on a different seed
# with more draws: 4 chains x 6000 (3000 warmup) on both sides, seed 101.
#   Rscript dev/sampleprior-rev-agree.R [case ...]
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sampleprior"
cache <- file.path(wt, "dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source(file.path(wt, "dev/sampleprior-agree-common.R"))
cat("lib:", find.package("frmtmb.sample"), "\n")
which <- commandArgs(TRUE)
cases <- agree_cases()
if (length(which)) cases <- cases[which]
chains <- 4L
iter <- 6000L
seed <- 101L
for (nm in names(cases)) {
  cs <- cases[[nm]]
  f <- agree_frm(cs, seed, chains, iter)
  b <- agree_brms(cs, seed, chains, iter, cache)
  lf <- log(f)
  lb <- log(b)
  mf <- agree_mcse(lf, chains)
  mb <- agree_mcse(lb, chains)
  ef <- posterior::ess_bulk(matrix(lf, ncol = chains))
  eb <- posterior::ess_bulk(matrix(lb, ncol = chains))
  z <- (mean(lf) - mean(lb)) / sqrt(mf^2 + mb^2)
  # quantiles too: the mean of log(nu) is sensitive to the upper tail
  qs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  qf <- stats::quantile(lf, qs)
  qb <- stats::quantile(lb, qs)
  mq <- vapply(qs, function(p) {
    c(posterior::mcse_quantile(matrix(lf, ncol = chains), probs = p),
      posterior::mcse_quantile(matrix(lb, ncol = chains), probs = p))
  }, c(0, 0))
  zq <- (qf - qb) / sqrt(colSums(mq^2))
  cat(sprintf(paste0("%-11s log(%s): frmtmb mean %.4f mcse %.4f ess %.0f | ",
                     "brms mean %.4f mcse %.4f ess %.0f | z(mean) %.2f | ",
                     "z(q05,q25,q50,q75,q95) %s | sd ratio %.3f\n"),
              nm, cs$par, mean(lf), mf, ef, mean(lb), mb, eb, z,
              paste(sprintf("%.2f", zq), collapse = " "),
              stats::sd(lf) / stats::sd(lb)))
}
