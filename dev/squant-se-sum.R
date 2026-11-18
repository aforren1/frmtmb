# Lane squant: summarize dev/squant-se.R's tables. For each design, n,
# quantile and method: the Monte Carlo sd of the slope and intercept,
# the mean Hessian and sandwich standard errors (rq: its "nid"
# sandwich), the bias against the true conditional quantile, and the
# coverage of 95% Wald intervals from each standard error.
#
# Run: Rscript dev/squant-se-sum.R > dev/squant-out/se-summary.txt
fs <- Sys.glob("dev/squant-out/se-*-*.tsv")
d <- do.call(rbind, lapply(fs, function(f) {
  x <- utils::read.delim(f, stringsAsFactors = FALSE, quote = "")
  x[, setdiff(names(x), "warnings")]
}))
w <- do.call(rbind, lapply(fs, function(f) {
  utils::read.delim(f, stringsAsFactors = FALSE, quote = "")
}))
d$h[is.na(d$h)] <- -1
w$h[is.na(w$h)] <- -1
key <- interaction(d$design, d$n, d$tau, d$method, d$h, drop = TRUE,
                   lex.order = TRUE)
cat("design\tn\ttau\tmethod\th\treps\tconv0\twarned\t",
    "bias1/sd\tsd1\tse_hess1/sd\tse_sand1/sd\tcov_hess1\tcov_sand1\t",
    "bias0/sd\tsd0\tse_hess0/sd\tse_sand0/sd\tcov_hess0\tcov_sand0\n",
    sep = "")
for (k in levels(key)) {
  x <- d[key == k & is.finite(d$b1), ]
  ww <- w[interaction(w$design, w$n, w$tau, w$method, w$h, drop = TRUE,
                       lex.order = TRUE) == k, ]
  if (!nrow(x)) next
  cov <- function(b, se, tr) mean(abs(b - tr) <= 1.959964 * se, na.rm = TRUE)
  sd1 <- sd(x$b1)
  sd0 <- sd(x$b0)
  cat(sprintf(paste0("%s\t%d\t%.2f\t%s\t%s\t%d\t%d\t%d\t%+.3f\t%.4f\t%.3f\t",
                     "%.3f\t%.3f\t%.3f\t%+.3f\t%.4f\t%.3f\t%.3f\t%.3f\t",
                     "%.3f\n"),
              x$design[1], x$n[1], x$tau[1], x$method[1],
              ifelse(x$h[1] < 0, "-", format(x$h[1])), nrow(x),
              sum(x$conv == 0), sum(nzchar(ww$warnings) & !is.na(ww$warnings)),
              (mean(x$b1) - x$true1[1]) / sd1, sd1,
              mean(x$se_hess1) / sd1, mean(x$se_sand1) / sd1,
              cov(x$b1, x$se_hess1, x$true1), cov(x$b1, x$se_sand1, x$true1),
              (mean(x$b0) - x$true0[1]) / sd0, sd0,
              mean(x$se_hess0) / sd0, mean(x$se_sand0) / sd0,
              cov(x$b0, x$se_hess0, x$true0), cov(x$b0, x$se_sand0, x$true0)))
}
wm <- w$warnings[!is.na(w$warnings) & nzchar(w$warnings)]
cat("\nwarning texts (count):\n")
print(sort(table(substr(wm, 1, 80)), decreasing = TRUE))
