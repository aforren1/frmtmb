# Lane `latent`, punch round 1: the lca recovery table under the start
# that now ships.
#
#   Rscript dev/latent-2p4-recheck-summarize.R [tsv]

args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args)) args[[1L]] else "dev/latent-2p4-recheck.tsv"
d <- utils::read.delim(f, check.names = FALSE)
R <- nrow(d)
cat("replicates:", R, "  file:", f, "\n")
cat("seeds:", min(d$seed), "..", max(d$seed), "\n\n")
q <- function(x) c(min = min(x), median = stats::median(x), max = max(x))

cat("== against poLCA(nrep = 10), read from the 200-replicate sweep ==\n")
cat("  reached poLCA's optimum on", sum(d$reached_polca), "of", R, "\n")
lo <- d$reached_polca == 0L
if (any(lo)) {
  cat("  the ones it did not:\n")
  print(d[lo, c("seed", "ll_frm", "ll_polca")], row.names = FALSE)
}
cat("  |ll_frm - ll_polca| / |ll_frm|  :",
    formatC(q(d$ll_rel_gap), digits = 3, format = "e"), "\n")

cat("\n== recovery against the simulator's truth ==\n")
cat("  item probabilities, max |error| :",
    formatC(q(d$prof_err_max), digits = 4, format = "g"), "\n")
cat("  the same table's RMSE           :",
    formatC(q(d$prof_err_rmse), digits = 4, format = "g"), "\n")
cat("  class shares, max |error|       :",
    formatC(q(d$share_err_max), digits = 4, format = "g"), "\n")
cat("  relative entropy                :",
    formatC(q(d$entropy), digits = 4, format = "g"), "\n")
cat("  modal-assignment accuracy       :",
    formatC(q(d$modal_acc), digits = 4, format = "g"), "\n")

nm <- c("c2_1", "c2_2", "c2_3", "c3_1", "c3_2", "c3_3",
        "c4_1", "c4_2", "c4_3")
truth <- c(-0.4, 0.8, -0.5, 0.2, -0.6, 0.9, -0.1, 0.3, 0.4)
lab <- c("class2:(Intercept)", "class2:x1", "class2:x2",
         "class3:(Intercept)", "class3:x1", "class3:x2",
         "class4:(Intercept)", "class4:x1", "class4:x2")
cat(sprintf("\n  %-20s %8s %9s %9s %9s %9s %8s\n", "coefficient",
            "truth", "bias", "RMSE", "emp.SD", "mean SE", "cover95"))
for (i in seq_along(nm)) {
  e <- d[[paste0("err_", nm[i])]]
  s <- d[[paste0("se_", nm[i])]]
  cat(sprintf("  %-20s %8.2f %9.4f %9.4f %9.4f %9.4f %7.1f%%\n",
              lab[i], truth[i], mean(e), sqrt(mean(e^2)), stats::sd(e),
              mean(s), 100 * mean(abs(e) < stats::qnorm(0.975) * s)))
}
p <- sum(d$coef_cover_9) / (9 * R)
se <- sqrt(p * (1 - p) / (9 * R))
cat(sprintf("\n  overall Wald coverage: %.1f%% (%d of %d),",
            100 * p, sum(d$coef_cover_9), 9 * R))
cat(sprintf(" Monte Carlo interval %.1f%% to %.1f%%\n",
            100 * (p - 1.96 * se), 100 * (p + 1.96 * se)))
cat("\n  positive definite Hessian in", sum(d$pdhess), "of", R, "\n")
cat("  max|grad| / |logLik|            :",
    formatC(q(d$maxgrad_rel), digits = 3, format = "e"), "\n")
cat("  frm() seconds                   :",
    formatC(q(d$s_frm), digits = 3, format = "g"), "\n")
