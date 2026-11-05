# Summary of the Laplace simulation. No model code: reads the sim tsv files.
args <- commandArgs(TRUE)
dir <- if (length(args)) args[[1]] else "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/sim"
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
truth <- list(
  pois = c(b0 = -0.2, b1 = 0.5, s = 0.5, ar = 0.6),
  pois_ar9 = c(b0 = -0.2, b1 = 0.5, s = 0.3, ar = 0.9),
  pois_s1 = c(b0 = -0.5, b1 = 0.5, s = 1, ar = 0.6),
  bern = c(b0 = 0.2, b1 = 0.5, s = 1, ar = 0.6),
  bern_long = c(b0 = 0.2, b1 = 0.5, s = 1, ar = 0.8),
  binom3 = c(b0 = 0.2, b1 = 0.5, s = 0.8, ar = 0.6),
  pois_big = c(b0 = 3, b1 = 0.3, s = 0.5, ar = 0.6))
cover <- function(est, se, tr, z = 1.96) mean(abs(est - tr) <= z * se, na.rm = TRUE)
for (dn in unique(x$design)) {
  d <- x[x$design == dn, ]
  tr <- truth[[dn]]
  cat("\n====", dn, " n =", nrow(d), " errors:", sum(nzchar(d$err) & !is.na(d$err)),
      " code1:", sum(d$code != 0, na.rm = TRUE), " warned:", sum(d$nwarn > 0, na.rm = TRUE),
      " code!=0 with no warning:", sum(d$code != 0 & d$nwarn == 0, na.rm = TRUE), "\n")
  mk <- function(lab, est, ex, t0) {
    data.frame(par = lab, truth = t0,
               frm_mean = mean(est, na.rm = TRUE), frm_median = median(est, na.rm = TRUE),
               exact_mean = mean(ex, na.rm = TRUE), exact_median = median(ex, na.rm = TRUE),
               med_diff_frm_exact = median(est - ex, na.rm = TRUE),
               frm_sd = sd(est, na.rm = TRUE), exact_sd = sd(ex, na.rm = TRUE))
  }
  tb <- rbind(mk("b0", d$b0, d$x_b0, tr["b0"]), mk("b1", d$b1, d$x_b1, tr["b1"]),
              mk("sderr", d$sderr, d$x_sderr, tr["s"]), mk("ar", d$ar, d$x_ar, tr["ar"]),
              mk("msd", d$sderr / sqrt(1 - pmin(d$ar^2, 1 - 1e-8)),
                 d$x_sderr / sqrt(1 - pmin(d$x_ar^2, 1 - 1e-8)),
                 tr["s"] / sqrt(1 - tr["ar"]^2)))
  print(tb, digits = 3, row.names = FALSE)
  th1 <- tr["ar"] / sqrt(1 - tr["ar"]^2)
  # coverage, Wald on the fit's reported scale
  cat(sprintf("coverage b0 %.3f  b1 %.3f  sderr(CI) %.3f  ar(CI) %.3f  log sderr(wald) %.3f  chart ar(wald) %.3f\n",
    cover(d$b0, d$se_b0, tr["b0"]), cover(d$b1, d$se_b1, tr["b1"]),
    mean(d$sderr_lo <= tr["s"] & d$sderr_hi >= tr["s"], na.rm = TRUE),
    mean(d$ar_lo <= tr["ar"] & d$ar_hi >= tr["ar"], na.rm = TRUE),
    cover(d$ls, d$se_ls, log(tr["s"])), cover(d$t1, d$se_t1, th1)))
  cat(sprintf("NA CI: sderr %d of %d, ar %d; NA SE: b0 %d, ls %d, t1 %d\n", sum(is.na(d$sderr_lo)), nrow(d), sum(is.na(d$ar_lo)), sum(is.na(d$se_b0)), sum(is.na(d$se_ls)), sum(is.na(d$se_t1))))
  # at the boundary
  cat(sprintf("frm stops with ar>0.99: %d; sderr<0.01: %d; exact stops at ar bound: %d, at s bound(>=9.9): %d, s<0.021: %d\n",
    sum(d$ar > 0.99, na.rm = TRUE), sum(d$sderr < 0.01, na.rm = TRUE),
    sum(abs(d$x_t1) >= 7.99, na.rm = TRUE), sum(d$x_sderr >= 9.9, na.rm = TRUE), sum(d$x_sderr < 0.021, na.rm = TRUE)))
  # Laplace error in the value, at frm()'s own point
  ok <- abs(d$t1) < 8 & !is.na(d$x_at_frm)
  cat(sprintf("Laplace - exact at frm's point: median %.4f, q05 %.4f, q95 %.4f (n=%d); exact max - exact at frm's point: median %.4f, q95 %.4f\n",
    median((d$ll - d$x_at_frm)[ok]), quantile((d$ll - d$x_at_frm)[ok], .05), quantile((d$ll - d$x_at_frm)[ok], .95), sum(ok),
    median((d$x_ll - d$x_at_frm)[ok]), quantile((d$x_ll - d$x_at_frm)[ok], .95)))
  cat("warn kinds:\n")
  w <- substr(d$warn[!is.na(d$warn) & nzchar(d$warn)], 1, 60)
  print(utils::head(sort(table(w), decreasing = TRUE), 5))
}
