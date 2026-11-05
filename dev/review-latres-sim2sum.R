# Summary of review-latres-sim2.R: frm() against the exact MLE and the truth.
dir <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/sim2"
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
truth <- list(gamma_T = c(b1 = .4, s = .5, ar = .5, shape = 4), nb_T = c(b1 = .4, s = .5, ar = .5, shape = 4),
              gamma_F = c(b1 = .4, s = .5, ar = .5, shape = 4), probit_F = c(b1 = .5, s = .5, ar = .2),
              pois_F = c(b1 = .4, s = .5, ar = .5), bern_F = c(b1 = .5, s = 1, ar = .5))
for (dn in names(truth)) {
  d <- x[x$design == dn, ]
  if (!nrow(d)) next
  tr <- truth[[dn]]
  cat("\n====", dn, " n", nrow(d), " errors", sum(nzchar(d$err) & !is.na(d$err)), " code!=0", sum(d$code != 0, na.rm = TRUE),
      " warned", sum(d$nwarn > 0, na.rm = TRUE), " code!=0 unwarned", sum(d$code != 0 & d$nwarn == 0, na.rm = TRUE),
      " mean tfit", round(mean(d$tfit), 2), " mean tex", round(mean(d$tex, na.rm = TRUE), 2), "\n")
  mk <- function(lab, est, ex, t0) data.frame(par = lab, truth = t0, frm_med = median(est, na.rm = TRUE), frm_mean = mean(est, na.rm = TRUE), exact_med = median(ex, na.rm = TRUE), exact_mean = mean(ex, na.rm = TRUE), med_diff = median(est - ex, na.rm = TRUE), frm_sd = sd(est, na.rm = TRUE), ex_sd = sd(ex, na.rm = TRUE))
  tb <- rbind(mk("b1", d$frm_b1, d$x_b1, tr["b1"]), mk("sderr", exp(d$frm_ls), exp(d$x_ls), tr["s"]), mk("ar", d$frm_ar, d$x_ar, tr["ar"]))
  if (!is.na(tr["shape"])) tb <- rbind(tb, mk("shape", exp(d$frm_lshape), exp(d$x_lshape), tr["shape"]))
  print(tb, digits = 3, row.names = FALSE)
  if (!is.na(tr["shape"])) {
    cat("shape > 1000: frm", sum(d$frm_lshape > log(1000), na.rm = TRUE), " exact", sum(d$x_lshape > log(1000), na.rm = TRUE), " both", sum(d$frm_lshape > log(1000) & d$x_lshape > log(1000), na.rm = TRUE), "of", nrow(d), "\n")
  }
  ok <- is.finite(d$ll) & is.finite(d$x_at_frm) & d$x_at_frm > -1e9
  cat(sprintf("logLik(frm) - exact at frm's point: median %.3f, q05 %.3f, q95 %.3f; exact max - frm logLik: median %.3f, q95 %.3f; exact max - exact at frm point: median %.3f, q95 %.3f, max %.3f (n=%d)\n",
    median((d$ll - d$x_at_frm)[ok]), quantile((d$ll - d$x_at_frm)[ok], .05), quantile((d$ll - d$x_at_frm)[ok], .95),
    median((d$x_ll - d$ll)[ok]), quantile((d$x_ll - d$ll)[ok], .95),
    median((d$x_ll - d$x_at_frm)[ok]), quantile((d$x_ll - d$x_at_frm)[ok], .95), max((d$x_ll - d$x_at_frm)[ok]), sum(ok)))
  cat("wald coverage b1:", mean(abs(d$frm_b1 - tr["b1"]) <= 1.96 * d$se_b1, na.rm = TRUE), " log sderr:", mean(abs(d$frm_ls - log(tr["s"])) <= 1.96 * d$se_ls, na.rm = TRUE), "\n")
  print(utils::head(sort(table(substr(d$warn[!is.na(d$warn) & nzchar(d$warn)], 1, 50)), decreasing = TRUE), 4))
}
