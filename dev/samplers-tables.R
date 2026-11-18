# Markdown tables for dev/samplers-findings.md, generated from the
# lane's per-design TSV rows.
#   Rscript dev/samplers-tables.R > dev/samplers-tables.md
rd <- function(dir) {
  f <- list.files(dir, "\\.tsv$", full.names = TRUE)
  do.call(rbind, lapply(f, function(x) read.delim(x, stringsAsFactors = FALSE)))
}
md <- function(d) {
  d[] <- lapply(d, function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    gsub("|", "\\|", x, fixed = TRUE)
  })
  cat("| ", paste(names(d), collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", ncol(d)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(d))) {
    cat("| ", paste(unlist(d[i, ]), collapse = " | "), " |\n", sep = "")
  }
  cat("\n")
}
ord <- function(d, first) d[order(!d$id %in% first, d$id), ]

q1 <- rd("dev/samplers-q1-out")
if (!is.null(q1)) {
  q1 <- q1[order(startsWith(q1$id, "big_"), q1$id), ]
  q1$stanli_total <- q1$t_stanli_setup + q1$t_stanli_sample
  q1$rstan_total <- q1$t_rstan_compile + q1$t_rstan_sample
  q1$speedup_sample <- round(q1$t_rstan_sample / q1$t_stanli_sample, 1)
  q1$speedup_total <- round(q1$rstan_total / q1$stanli_total, 1)
  cat("### Q1a. Wall time (s), 4 chains x 1000/1000 run one after another\n\n")
  md(q1[, c("id", "key", "n", "t_frm", "t_codegen", "t_stanli_setup",
            "t_stanli_sample", "t_rstan_compile", "t_rstan_sample",
            "speedup_sample", "speedup_total")])
  cat("### Q1b. stanli against rstan: diagnostics and posterior means\n\n")
  cat("z = (mean_stanli - mean_rstan) / sqrt(mcse_stanli^2 + mcse_rstan^2)",
      "over every variable both report (r_ included, lp__ included);",
      "crit = the two-sided 5% Bonferroni bound for that many z.\n\n")
  md(q1[, c("id", "stanli_div", "rstan_div", "stanli_maxrhat",
            "rstan_maxrhat", "stanli_miness", "rstan_miness", "sr_nvar",
            "sr_maxz", "sr_crit", "sr_frac3", "sr_lp_z", "stanli_error",
            "rstan_error")])
  cat("### Q1c. frm() maximum likelihood against the stanli posterior",
      "(brms default priors)\n\n")
  cat("z = (MLE - posterior median) / posterior sd. `mle_maxz_hyp` over",
      "b_, sd_, cor_, sigma and other population parameters;",
      "`mle_medz_re` = median |z| over r_ levels.\n\n")
  md(q1[, c("id", "mle_nhyp", "mle_maxz_hyp", "mle_worst", "mle_medz_re",
            "frm_error", "mle_error")])
}

q2 <- rd("dev/samplers-q2-out")
if (!is.null(q2)) {
  cat("### Q2a. Survey instrument steps (s) and the no-Jacobian shim\n\n")
  cat("Every program was in dev/stan-cache, so `t_rstan_compile` is",
      "rstan reusing the program it had just read, not a compile; the",
      "cold compile times are Q1a's.",
      "`d_*` = max relative difference, stanli shim against rstan, at",
      "frm()'s estimate and 5 points around it.\n\n")
  md(q2[, c("id", "key", "npar_u", "cache_hit", "t_rstan_warm",
            "t_rstan_chains0", "t_rstan_lp",
            "t_stanli_setup", "t_stanli_setup_jac", "t_stanli_lp", "d_noj",
            "d_j", "d_grad", "d_unconstrain")])
  cat("### Q2b. L4 sampling (2 chains x 2000, flat priors) and stanli's",
      "rstan-shaped calls\n\n")
  md(q2[, c("id", "t_rstan_l4", "t_stanli_l4", "l4_nvar", "l4_maxz",
            "l4_crit", "l4_rhat_rstan", "l4_rhat_stanli", "asm_chains0",
            "asm_noj", "fp_error")])
}
fp <- "dev/samplers-q2-out/fp.tsv.txt"
if (file.exists(fp)) {
  cat("### Q2c. The method tier's Fixed_param brmsfit on as_stanmodel()",
      "(keep-list with include = TRUE)\n\n")
  md(read.delim(fp, stringsAsFactors = FALSE))
}

q3 <- rd("dev/samplers-q3-out")
if (!is.null(q3)) {
  q3 <- ord(q3, c("small_pois", "big_pois", "big_crossed", "slope", "gp"))
  # min bulk ESS over the variables BOTH arms report: the tmbstan draws
  # carry lp__, the SNUTS draws as mapped do not, so the TSV's minimum
  # could be lp__ on one side only
  for (i in seq_len(nrow(q3))) {
    r <- readRDS(file.path("dev/samplers-q3-out", paste0(q3$id[i], ".rds")))
    a <- r$sum_tmb
    b <- r$sum_sn
    if (is.null(a) || is.null(b)) next
    v <- setdiff(intersect(a$variable, b$variable), "lp__")
    q3$tmb_miness[i] <- round(min(a$ess_bulk[match(v, a$variable)]), 1)
    q3$sn_miness[i] <- round(min(b$ess_bulk[match(v, b$variable)]), 1)
    q3$tmb_ess_s[i] <- round(q3$tmb_miness[i] / r$t_tmb_total, 3)
    q3$sn_ess_s[i] <- round(q3$sn_miness[i] / r$t_sn_total, 3)
  }
  q3$ratio <- round(q3$sn_ess_s / q3$tmb_ess_s, 2)
  cat("### Q3a. Efficiency: min bulk ESS over every variable per second",
      "of wall time (whole call)\n\n")
  md(q3[, c("id", "n", "n_par", "n_random", "tmb_total", "tmb_sampling",
            "tmb_miness", "tmb_ess_s", "sn_metric", "sn_total",
            "sn_sampling", "sn_miness", "sn_ess_s", "ratio")])
  cat("### Q3b. Diagnostics and agreement (tmbstan route vs SNUTS)\n\n")
  cat("z as in Q1b over every variable both arms report;",
      "sd ratio = posterior sd tmbstan / SNUTS.\n\n")
  md(q3[, c("id", "tmb_div", "tmb_maxtd", "tmb_maxrhat", "tmb_leap",
            "sn_div", "sn_maxtd", "sn_maxrhat", "sn_leap", "nvar", "maxz",
            "maxz_hyp", "crit", "worst", "sdratio_range", "error")])
}
