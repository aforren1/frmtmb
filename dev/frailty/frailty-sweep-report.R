source("frailty-lib.R")
x <- readLines("frailty-sweep.tsv")
rows <- lapply(strsplit(x, "\t"), function(f) {
  kv <- strsplit(f, "=")
  v <- vapply(kv, function(z) trimws(z[2L]), character(1))
  names(v) <- vapply(kv, `[`, character(1), 1L)
  v
})
nm <- names(rows[[1L]])
d <- as.data.frame(do.call(rbind, lapply(rows, function(r) r[nm])),
                   stringsAsFactors = FALSE)
for (j in setdiff(nm, "cell")) d[[j]] <- as.numeric(d[[j]])
cells <- unique(d$cell)
f <- function(v) formatC(v, digits = 3, format = "g", width = 10)
cat(sprintf("%-8s %6s %6s %6s %10s %10s %10s %10s %10s %10s %10s\n",
            "cell", "ncent", "per", "ev/c", "beta_rel", "beta/se",
            "sd_rel", "lapErr", "ghqErr", "exactGap", "map_rel"))
for (cl in cells) {
  s <- d[d$cell == cl, ]
  cat(sprintf("%-8s %6d %6d %6.1f %s %s %s %s %s %s %s\n", cl,
              s$n_centre[1L], s$per_centre[1L], mean(s$ev_per_centre),
              f(mean(s$beta_reldiff)), f(mean(s$beta_diff_in_se)),
              f(mean(s$sd_reldiff)), f(mean(s$lap_err_at_frm)),
              f(mean(s$ghq_err_at_rst)), f(mean(s$exact_gap)),
              f(max(s$map_resid_rel))))
}
cat("\nworst single replicate per cell (max over 5):\n")
cat(sprintf("%-8s %10s %10s %10s %10s\n", "cell", "beta/se", "sd_rel",
            "lapErr", "exactGap"))
for (cl in cells) {
  s <- d[d$cell == cl, ]
  cat(sprintf("%-8s %s %s %s %s\n", cl, f(max(s$beta_diff_in_se)),
              f(max(s$sd_reldiff)), f(min(s$lap_err_at_frm)),
              f(min(s$exact_gap))))
}
cat("\nse agreement and coverage (5 replicates per cell):\n")
cat(sprintf("%-8s %12s %12s %12s %12s %8s %8s\n", "cell",
            "se_b_frm", "se_b_rst", "se_lsd_frm", "se_lsd_rst",
            "b_cov", "sd_cov"))
for (cl in cells) {
  s <- d[d$cell == cl, ]
  cat(sprintf("%-8s %12s %12s %12s %12s %8s %8s\n", cl,
              f(mean(s$se_beta_frm)), f(mean(s$se_beta_rst)),
              f(mean(s$se_logsd_frm)), f(mean(s$se_logsd_rst)),
              paste0(sum(s$b_cover), "/5"), paste0(sum(s$sd_cover), "/5")))
}
cat("\ntiming, seconds per fit (mean of 5):\n")
for (cl in cells) {
  s <- d[d$cell == cl, ]
  cat(sprintf("%-8s frmtmb %6.2f  rstpm2 %6.2f\n", cl, mean(s$t_frm),
              mean(s$t_rst)))
}
