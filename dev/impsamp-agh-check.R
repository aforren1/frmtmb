.libPaths(c("C:/Users/adf44/source/r/rellib-r16","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("dev/impsamp-agh.R")
for (fam in c("binomial","poisson","negbinomial","gaussian","Gamma","bernoulli")) {
  for (sh in c("int","slope")) {
  d <- impsamp_data(fam, sh, "10x20", 1)
  form <- impsamp_formula(fam, sh)
  lap <- suppressWarnings(suppressMessages(frm(form, family=impsamp_family(fam), data=d)))
  st <- impsamp_start(fam, sh, lap)
  cat(fam, sh, "par", format(lap$opt$par, digits=5), "\n start", format(st,digits=5), "\n")
  t0 <- proc.time()[["elapsed"]]
  e <- agh_ll(fam, sh, d, st, if (sh=="int") 25L else 11L)
  el <- proc.time()[["elapsed"]]-t0
  q <- if (sh=="int") suppressWarnings(suppressMessages(frm(form, family=impsamp_family(fam), data=d, quadrature=TRUE)))
  cat(sprintf(" lap %.5f  agh@lap %.5f (%.2fs)  quad %s quad's agh@own %s\n", as.numeric(logLik(lap)), e, el,
    if (is.null(q)) "NA" else format(as.numeric(logLik(q)), digits=10),
    if (is.null(q)) "NA" else format(agh_ll(fam, sh, d, impsamp_start(fam, sh, q), 25L), digits=10)))
  }
}
