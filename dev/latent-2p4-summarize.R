# Lane `latent`, item 2.4: the recovery table and the identity table
# from dev/latent-2p4-lca.tsv.
#
#   Rscript dev/latent-2p4-summarize.R [tsv]

args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args)) args[[1L]] else "dev/latent-2p4-lca.tsv"
d <- utils::read.delim(f, check.names = FALSE)
R <- nrow(d)
cat("replicates:", R, "  file:", f, "\n")
cat("seeds:", min(d$seed), "..", max(d$seed),
    "  poLCA start seeds:", min(d$seed_polca), "..", max(d$seed_polca),
    "\n\n")

q <- function(x) {
  c(min = min(x), median = stats::median(x), max = max(x))
}

cat("== identity against poLCA, same data, both optimized ==\n")
cat("  |ll_frm - ll_poLCA| / |ll_frm|      :",
    formatC(q(d$ll_rel_gap), digits = 3, format = "e"), "\n")
cat("  replicates where frm reached a HIGHER logLik than poLCA's",
    "best of 10 EM starts:", sum(d$frm_beats_polca), "of", R, "\n")
cat("  worst signed gap ll_frm - ll_poLCA  :",
    formatC(max(abs(d$ll_rel_gap)) * max(abs(as.numeric(d$ll_frm))),
            digits = 3, format = "e"), "\n")
cat("  item profiles, max |difference|     :",
    formatC(q(d$prof_id_max), digits = 3, format = "e"), "\n")
cat("    as a ratio to poLCA's own largest profile SE:",
    formatC(q(d$prof_id_over_polca_se), digits = 3, format = "e"), "\n")
cat("  posterior class probs, max |diff|   :",
    formatC(q(d$post_id_max), digits = 3, format = "e"), "\n")
cat("  gating coefficients, max |diff|     :",
    formatC(q(d$coef_id_max), digits = 3, format = "e"), "\n")
cat("    as a ratio to poLCA's own largest coefficient SE:",
    formatC(q(d$coef_id_over_polca_se), digits = 3, format = "e"), "\n")

cat("\n== recovery against the simulator's truth ==\n")
cat("  item endorsement probabilities, max |error| over the 4 x 10",
    "table:\n     ", formatC(q(d$prof_err_max), digits = 4,
                             format = "g"), "\n")
cat("  the same table's RMSE               :",
    formatC(q(d$prof_err_rmse), digits = 4, format = "g"), "\n")
cat("  class shares, max |error|           :",
    formatC(q(d$share_err_max), digits = 4, format = "g"), "\n")
cat("  relative entropy of the classification:",
    formatC(q(d$entropy), digits = 4, format = "g"), "\n")
cat("  modal-assignment accuracy           :",
    formatC(q(d$modal_acc), digits = 4, format = "g"), "\n")

nm <- c("c2_1", "c2_2", "c2_3", "c3_1", "c3_2", "c3_3",
        "c4_1", "c4_2", "c4_3")
truth <- c(-0.4, 0.8, -0.5, 0.2, -0.6, 0.9, -0.1, 0.3, 0.4)
lab <- c("class2:(Intercept)", "class2:x1", "class2:x2",
         "class3:(Intercept)", "class3:x1", "class3:x2",
         "class4:(Intercept)", "class4:x1", "class4:x2")
cat("\n  gating coefficients, re-referenced to class 1 and matched to",
    "the truth by profile:\n")
cat(sprintf("  %-20s %8s %9s %9s %9s %9s %9s\n", "coefficient", "truth",
            "bias", "RMSE", "emp.SD", "mean SE", "cover95"))
for (i in seq_along(nm)) {
  e <- d[[paste0("err_", nm[i])]]
  s <- d[[paste0("se_", nm[i])]]
  cov <- mean(abs(e) < stats::qnorm(0.975) * s)
  cat(sprintf("  %-20s %8.2f %9.4f %9.4f %9.4f %9.4f %8.1f%%\n",
              lab[i], truth[i], mean(e), sqrt(mean(e^2)), stats::sd(e),
              mean(s), 100 * cov))
}
cat(sprintf("\n  overall Wald coverage over %d replicates x 9",
            R))
cat(sprintf(" coefficients: %.1f%% (%d of %d)\n",
            100 * sum(d$coef_cover_9) / (9 * R), sum(d$coef_cover_9),
            9 * R))
# a binomial interval on the coverage itself, so the number is read with
# its own uncertainty rather than as an exact rate
p <- sum(d$coef_cover_9) / (9 * R)
se <- sqrt(p * (1 - p) / (9 * R))
cat(sprintf("  (its own Monte Carlo 95%% interval, treating the 9 x %d",
            R))
cat(sprintf(" as independent: %.1f%% to %.1f%%)\n",
            100 * (p - 1.96 * se), 100 * (p + 1.96 * se)))

cat("\n== diagnostics and cost ==\n")
cat("  positive definite Hessian in", sum(d$pdhess), "of", R, "\n")
cat("  max|grad| / |logLik|                :",
    formatC(q(d$maxgrad_rel), digits = 3, format = "e"), "\n")
cat("  frm() seconds                       :",
    formatC(q(d$s_frm), digits = 3, format = "g"), "\n")
cat("  poLCA seconds (nrep = 10)           :",
    formatC(q(d$s_polca), digits = 3, format = "g"), "\n")
cat("  class permutations used (frm label -> true label), tally:\n")
print(table(d$perm_truth))
