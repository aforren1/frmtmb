# Lane `latent`, item 2.3: the recovery table and the identity table
# from dev/latent-2p3-hmm-A.tsv or -B.tsv.
#
#   Rscript dev/latent-2p3-summarize.R <tsv> [reps]
#
# `reps` truncates to the first N replicates, which is how a declared
# replicate count is honoured when the run was launched with a larger
# one: the first N are the first N SEEDS, so nothing is selected on its
# outcome.

args <- commandArgs(trailingOnly = TRUE)
f <- args[[1L]]
d <- utils::read.delim(f, check.names = FALSE)
if (length(args) >= 2L) {
  n <- min(nrow(d), as.integer(args[[2L]]))
  d <- d[seq_len(n), , drop = FALSE]
}
R <- nrow(d)
arm <- d$arm[1L]
cat("file:", f, "  arm:", arm, "  replicates:", R, "\n")
cat("seeds:", min(d$seed), "..", max(d$seed), "  rows per replicate:",
    unique(d$rows), "\n\n")

q <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(min = NA_real_, median = NA_real_, max = NA_real_))
  }
  c(min = min(x), median = stats::median(x), max = max(x))
}
fq <- function(x, dg = 3) formatC(q(x), digits = dg, format = "e")
gq <- function(x, dg = 4) formatC(q(x), digits = dg, format = "g")

ref <- if (identical(arm, "A")) "depmixS4 (best of 4 random EM starts)" else
  "hmmTMB (started at the truth)"
gapv <- suppressWarnings(as.numeric(d$ll_rel_gap))
nid <- sum(is.finite(gapv))
if (!nid) {
  cat("== no third-party arm in this file: recovery only ==\n")
} else {
  cat("== identity against", ref, "on the first", nid,
      "replicates ==\n")
  cat("  |ll_frm - ll_ref| / |ll_frm|   :", fq(gapv), "\n")
  big <- is.finite(gapv) & gapv > 1e-8
  cat("  replicates disagreeing by more than 1e-8 relative:", sum(big),
      "of", nid, "\n")
  if (any(big)) {
    cat("    their seeds:", d$seed[big], "\n")
    cat("    their signed gaps (frm minus reference):",
        formatC(d$ll_frm[big] - d$ll_ref[big], digits = 5,
                format = "f"), "\n")
  }
  cat("  max |estimate difference| (means, sds, sd(re)):",
      fq(d$est_max_gap), "\n")
  cat("  max |transition matrix difference|            :",
      fq(d$tpm_max_gap), "\n")
}

cat("\n== recovery against the simulator's truth ==\n")
cat("  state means, max |error|      :", gq(d$mu_err_max), "\n")
cat("  state sds,   max |error|      :", gq(d$sd_err_max), "\n")
cat("  transition matrix, max |error|:", gq(d$tpm_err_max), "\n")
if (identical(arm, "B")) {
  cat("  sd(tr12 | id), estimated      :", gq(d$sd_tr12),
      "   truth 0.6\n")
  cat("    mean", formatC(mean(d$sd_tr12), digits = 4, format = "g"),
      " sd", formatC(stats::sd(d$sd_tr12), digits = 4, format = "g"),
      "\n")
}
cat("  fitted labels were the true labels in", sum(d$identity_ok),
    "of", R, "replicates\n")
if (sum(d$identity_ok) < R) {
  cat("    permutations seen:", paste(unique(d$perm), collapse = " "),
      "\n")
}

nm <- c(paste0("mu", 1:3), paste0("sigma", 1:3),
        "tr12", "tr13", "tr22", "tr23", "tr32", "tr33",
        if (identical(arm, "B")) "theta1")
truth <- c(-2, 0, 3, rep(log(0.7), 3), -1.5, -2.5, 2.0, -1.0, -1.0, 2.0,
           if (identical(arm, "B")) log(0.6))
cat("\n  every parameter on its own internal scale, over the",
    sum(d$identity_ok), "replicates whose labels matched:\n")
cat(sprintf("  %-8s %8s %9s %9s %9s %9s %8s\n", "par", "truth", "bias",
            "RMSE", "emp.SD", "mean SE", "cover95"))
tot_c <- 0L
tot_n <- 0L
for (i in seq_along(nm)) {
  e <- d[[paste0("err_", nm[i])]]
  s <- d[[paste0("se_", nm[i])]]
  ok <- is.finite(e) & is.finite(s)
  if (!any(ok)) next
  cv <- d[[paste0("cov_", nm[i])]][ok]
  tot_c <- tot_c + sum(cv)
  tot_n <- tot_n + length(cv)
  cat(sprintf("  %-8s %8.3f %9.4f %9.4f %9.4f %9.4f %7.1f%%\n",
              nm[i], truth[i], mean(e[ok]), sqrt(mean(e[ok]^2)),
              stats::sd(e[ok]), mean(s[ok]), 100 * mean(cv)))
}
p <- tot_c / tot_n
se <- sqrt(p * (1 - p) / tot_n)
cat(sprintf("\n  overall Wald coverage: %.1f%% (%d of %d),",
            100 * p, tot_c, tot_n))
cat(sprintf(" Monte Carlo interval %.1f%% to %.1f%%\n",
            100 * (p - 1.96 * se), 100 * (p + 1.96 * se)))

cat("\n== diagnostics and cost ==\n")
cat("  positive definite Hessian in", sum(d$pdhess), "of", R, "\n")
cat("  max|grad| / |logLik|          :", fq(d$maxgrad_rel), "\n")
cat("  nlminb gradient calls         :", gq(d$nlminb_gr), "\n")
cat("  frm() seconds  (INDICATIVE, other work shared the machine):",
    gq(d$s_frm), "\n")
cat("  reference seconds (same caveat)                           :",
    gq(d$s_ref), "\n")
