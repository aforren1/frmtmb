# Reviewer of lane ordmix, re-check: the hurdle mixture with a prior on
# hu1 alone (dev/ordmix-rev2-b1.R case hu1_only), seeds 2, 7 and 8, whose
# full covariance has a standard error near 3e3 with no warning. Which
# parameter, and what the degenerate check measured.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ex <- parse("C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-rev-falsealarm.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("cut4", "gen")) eval(e)
}
for (s in c(2, 7, 8)) {
  d <- gen(s)
  f <- frm(bf(yh ~ x), family = mixture(hurdle_cumulative(),
                                        hurdle_cumulative()),
           data = d, prior = set_prior("beta(4, 16)", class = "hu1"))
  v <- sqrt(diag(vcov(f, full = TRUE)))
  cat("seed", s, "largest SEs:\n")
  print(round(head(sort(v, decreasing = TRUE), 3), 2))
  print(round(f$opt$par[names(v)[order(-v)][1:3]], 2))
  print(frmtmb:::mixture_ord_degeneracy(f, names(f$spec$responses)[1]))
}
