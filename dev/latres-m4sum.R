# Summary of dev/latres-m4.R: per design, frm() against the exact marginal, and
# the objective's dependence on the inner start.
fs <- list.files("dev/latres-out/m4", "[.]tsv$", full.names = TRUE)
o <- do.call(rbind, lapply(fs, function(f) {
  x <- read.delim(f, stringsAsFactors = FALSE)
  if (!nrow(x)) return(NULL)
  ok <- is.finite(x$ll)
  gap <- x$ll - x$x_at_frm
  data.frame(design = x$design[1], n = nrow(x), fits = sum(ok),
             code1 = sum(x$code[ok] != 0, na.rm = TRUE),
             silent = sum(x$code[ok] == 0 & x$nwarn[ok] == 0, na.rm = TRUE),
             sderr = median(x$sderr, na.rm = TRUE), x_sderr = median(x$x_sderr, na.rm = TRUE),
             ar = median(x$ar, na.rm = TRUE), x_ar = median(x$x_ar, na.rm = TRUE),
             gap_med = median(gap, na.rm = TRUE),
             gap_max = max(abs(gap), na.rm = TRUE),
             best_gap = median(x$x_ll - x$ll, na.rm = TRUE),
             span_med = median(x$span, na.rm = TRUE), span_max = max(x$span, na.rm = TRUE),
             nan_starts = median(x$nan_starts, na.rm = TRUE),
             ll_max = max(x$ll, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
print(format(o, digits = 3), row.names = FALSE)
