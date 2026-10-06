# Punch round 2, B3: the two suite fits the round-2 rule leaves silent
# (dev/ordmix-p2-log-suitefires.txt), both with cs(): what the
# component statistics are, and the estimates.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- parse("C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-p1-suitefires.R")
for (e in src) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("omx_data", "alt_data")) eval(e)
}
fits <- list(
  refusals_cs = frm(bf(y ~ x, mu2 ~ cs(x)),
                    family = mixture(cumulative(), sratio()),
                    data = omx_data(20261012)),
  brms_cs = frm(bf(y ~ cs(x)), family = mixture(sratio(), acat()),
                data = alt_data(20261030, 300)))
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("==", nm, "\n")
  print(frmtmb:::mixture_ord_degeneracy(f, names(f$spec$responses)[1]))
  print(suppressWarnings(fixef(f)))
}
