# Summary of dev/latres-m1.R: per design, frm() against the exact MLE, and the
# statistics of the fit a warning could read.
fs <- list.files("dev/latres-out/m1", "[.]tsv$", full.names = TRUE)
o <- do.call(rbind, lapply(fs, function(f) {
  x <- read.delim(f, stringsAsFactors = FALSE)
  if (!nrow(x)) return(NULL)
  bound <- mean(x$x_sderr >= 9.9, na.rm = TRUE)
  r <- x$sderr / pmin(x$x_sderr, 10)
  data.frame(design = x$design[1], n = nrow(x), rows = x$n[1],
             sderr = median(x$sderr, na.rm = TRUE),
             x_sderr = median(x$x_sderr, na.rm = TRUE),
             ratio = median(r, na.rm = TRUE),
             ar = median(x$ar, na.rm = TRUE), x_ar = median(x$x_ar, na.rm = TRUE),
             b1 = median(x$b1, na.rm = TRUE), x_b1 = median(x$x_b1, na.rm = TRUE),
             r_b1 = median(x$b1 / x$x_b1, na.rm = TRUE),
             at_bound = bound,
             lapl_gap = median(x$x_ll - x$ll, na.rm = TRUE),
             code1 = mean(x$code != 0, na.rm = TRUE),
             info_med = median(x$info_med, na.rm = TRUE),
             info_min = min(x$info_med, na.rm = TRUE),
             info_max = max(x$info_med, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
o <- o[order(o$info_med), ]
print(format(o, digits = 3), row.names = FALSE)
