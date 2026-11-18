# Short check of lane nlband: the kept hole fit of student() on gaussian
# band data, seed 2 (nu at its limit at the maximum): its standard
# errors and conditions, against the direct fit from (1, -2).
# NLBAND_ARM=after Rscript dev/reviews/nlband-short-studg.R x
source("dev/nlband-common.R")
d <- band_data(2)
r <- tt(frm(band_form, data = d, family = student()))
g <- tt(frm(band_form, data = d, family = student(),
            start = list(beta = c(1, -2)),
            control = frmtmb_control(nl_starts = FALSE)))
for (v in list(list("default", r), list("from (1,-2)", g))) {
  f <- v[[2]]
  se <- tt(sqrt(diag(vcov(f))))
  cat(sprintf("STUDG2 %s | par %s | se %s | fit msg %s | vcov() warn %s msg %s\n",
              v[[1]], paste(sprintf("%.5g", f$opt$par), collapse = " "),
              paste(signif(as.numeric(se), 4), collapse = " "),
              paste(attr(f, "mc"), collapse = "+"),
              paste(attr(se, "wc"), collapse = "+"),
              paste(attr(se, "mc"), collapse = "+")))
  for (m in attr(f, "m")) cat("   MSG", substr(gsub("\n", " ", m), 1, 160), "\n")
}
cat("DONE\n")
