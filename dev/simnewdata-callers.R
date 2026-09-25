# What the two callers that pass re_formula to simulate() do after the
# change: dharma_residuals() (default NULL) and pp_check() on a fit
# (default NA), on a mixed model and on a smooth, seeded. Run on both
# builds and compare the saved values.
#   SIMNEWDATA_LIB=base Rscript dev/simnewdata-callers.R base
#   Rscript dev/simnewdata-callers.R lane
#   Rscript dev/simnewdata-callers.R compare
arm <- commandArgs(trailingOnly = TRUE)[1]
OUT <- "dev/simnewdata-log"
if (identical(arm, "compare")) {
  a <- readRDS(file.path(OUT, "callers-base.rds"))
  b <- readRDS(file.path(OUT, "callers-lane.rds"))
  for (k in names(a)) {
    cat(sprintf("%-34s %s\n", k,
                if (identical(a[[k]], b[[k]])) "identical" else "differs"))
  }
  quit(save = "no")
}
source("dev/simnewdata-prelude.R")
suppressMessages(library(frmtmb))
set.seed(3)
d <- data.frame(x = stats::rnorm(200), w = stats::runif(200),
                g = factor(rep(1:10, 20)))
d$y <- 1 + 0.5 * d$x + stats::rnorm(10)[d$g] + stats::rnorm(200, 0, 0.5)
d$ys <- 2 * sin(2 * pi * d$w) + stats::rnorm(200, 0, 0.3)
fm <- frm(bf(y ~ x + (1 | g)), data = d)
fs <- frm(bf(ys ~ s(w)), data = d)
res <- list()
dh <- function(f, ...) {
  r <- dharma_residuals(f, nsim = 100, seed = 1, ...)
  r$scaledResiduals
}
res[["dharma mixed default"]] <- dh(fm)
res[["dharma mixed NA"]] <- dh(fm, re_formula = NA)
res[["dharma mixed ~1"]] <- dh(fm, re_formula = ~1)
res[["dharma smooth NA"]] <- dh(fs, re_formula = NA)
ppd <- function(f, ...) {
  set.seed(2)
  pp_check(f, ndraws = 20, ...)$data
}
res[["pp_check mixed default"]] <- ppd(fm)
res[["pp_check mixed NULL"]] <- ppd(fm, re_formula = NULL)
res[["pp_check mixed ~1"]] <- ppd(fm, re_formula = ~1)
res[["pp_check smooth default"]] <- ppd(fs)
# the smooth check's draws against the data: the spread of the
# replicated values around the observed ones
p <- ppd(fs)
yr <- p$value[!p$is_y]
cat(sprintf("%s: smooth pp_check, sd of replicated values %.3f, of y %.3f\n",
            arm, stats::sd(yr), stats::sd(d$ys)))
u <- res[["dharma smooth NA"]]
cat(sprintf("%s: smooth dharma(re_formula = NA), KS p against uniform %.3g\n",
            arm, stats::ks.test(u, "punif")$p.value))
u <- res[["dharma mixed ~1"]]
cat(sprintf("%s: mixed dharma(re_formula = ~1), KS p %.3g\n",
            arm, stats::ks.test(u, "punif")$p.value))
saveRDS(res, file.path(OUT, paste0("callers-", arm, ".rds")))
