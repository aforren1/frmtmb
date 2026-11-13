# Short check of lane ordstart's Q4 change: the five Q4 replacements of
# the sweep. Per fit: the warnings frm() gives, the thresholds of each
# component, the smallest probability any row gives each category per
# component, and a refit from the returned point (starts = 0): does the
# optimizer climb from there?
#   Rscript ordstart-short-q4pts.R
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/reviews/ordstart-review-sweep.R")
eval(parse(text = src[grep("^gen <- ", src):(grep("^cases <- list", src) - 1L)]))
q <- function(e) suppressWarnings(suppressMessages(e))
cases <- list(
  cs_sr_acat_seed4 = list(seed = 4, f = bf(y ~ cs(x) + z),
                          fam = mixture(sratio(), acat())),
  cs_sr_acat_seed1 = list(seed = 1, f = bf(y ~ cs(x) + z),
                          fam = mixture(sratio(), acat())),
  cs_sr_acat_seed3 = list(seed = 3, f = bf(y ~ cs(x) + z),
                          fam = mixture(sratio(), acat())),
  dpar_formula_seed4 = list(seed = 4, f = bf(y ~ 1, mu1 ~ x, mu2 ~ x + z),
                            fam = mixture(cumulative(), sratio())),
  dpar_formula_seed7 = list(seed = 7, f = bf(y ~ 1, mu1 ~ x, mu2 ~ x + z),
                            fam = mixture(cumulative(), sratio())))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  d <- gen(cs$seed, n = 300, K = 4)
  w <- character()
  f <- withCallingHandlers(frm(cs$f, family = cs$fam, data = d),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  cat(sprintf("\n=== %s ll %.4f code %d\n", nm, as.numeric(logLik(f)),
              f$opt$convergence))
  for (x in w) cat("  WARN:", substr(x, 1, 260), "\n")
  e <- f$estimates
  for (k in 1:2) {
    tr <- e[[paste0("tau_raw", k)]]
    cat(sprintf("  tau_raw%d: %s\n", k, paste(signif(tr, 5), collapse = " ")))
  }
  fe <- fixef(f)
  big <- fe[order(-abs(fe[, "Estimate"]))[1:6], c("Estimate", "Est.Error"),
            drop = FALSE]
  print(signif(big, 4))
  g <- q(frm(cs$f, family = cs$fam, data = d, start = f$estimates,
             control = frmtmb_control(starts = 0)))
  cat(sprintf("  REFIT from the returned point: ll %.4f code %d (gain %.4f)\n",
              as.numeric(logLik(g)), g$opt$convergence,
              as.numeric(logLik(g)) - as.numeric(logLik(f))))
}
cat("DONE\n")
