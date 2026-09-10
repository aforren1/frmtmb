# Reviewer, punch round 2, item 1: the under-coverage on the fresh
# block, recomputed from the raw columns, and asked of the TUNING block
# as well.
#
# The lane reports 93.4 percent at 1682 of 1800 on the fresh block
# against 94.3 percent on the tuning block, and attributes it to two of
# the nine gating coefficients. Two things the lane did not answer:
#
#   Q1  are the SAME two under-covering on the TUNING block? If they
#       are, this is a property of the design and not of the fresh
#       draw, and the help page should say so.
#   Q2  which coefficients are they, exactly, and how far under is the
#       reported standard error on each block?
#
# Nothing is fitted here. Everything comes from the lane's own
# permutation-matched `err_` and `se_` columns, which is the path this
# review used on the first two tables.
#
#   Rscript dev/rev-latent-cov3.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env3.R")
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
rd <- function(f) utils::read.delim(paste0(D, f), stringsAsFactors = FALSE)
num <- function(v) suppressWarnings(as.numeric(trimws(as.character(v))))
hr <- function(s) cat("\n==== ", s, " ====\n", sep = "")
mci <- function(k, n) {
  p <- k / n
  se <- sqrt(p * (1 - p) / n)
  sprintf("%.4f (%d/%d)  MC 95%%: %.4f to %.4f", p, k, n,
          p - 1.96 * se, p + 1.96 * se)
}

truth <- c(c2_1 = -0.4, c2_2 = 0.8, c2_3 = -0.5,
           c3_1 = 0.2, c3_2 = -0.6, c3_3 = 0.9,
           c4_1 = -0.1, c4_2 = 0.3, c4_3 = 0.4)
# the design's three gating columns, in the order the harness writes
# them: (Intercept), x1 (continuous), x2 (the binary covariate)
term <- c("(Intercept)", "x1", "x2 (binary)")

per_block <- function(f, label) {
  r <- rd(f)
  hr(paste(label, "-", f))
  cat("rows:", nrow(r), " seeds:", min(r$seed), "to", max(r$seed),
      " contiguous:",
      identical(sort(r$seed), min(r$seed):max(r$seed)), "\n")
  cat("reached poLCA (|d|/|ll| < 1e-6):",
      sum(abs(num(r$ll_rel_gap)) < 1e-6, na.rm = TRUE), "of", nrow(r),
      "  worst |d|/|ll|:",
      format(max(abs(num(r$ll_rel_gap)), na.rm = TRUE), digits = 4), "\n")
  cat("pdHess TRUE:", sum(num(r$pdhess) == 1 | r$pdhess %in% c(TRUE, "TRUE")),
      "of", nrow(r), "\n")
  cat("median fit seconds:",
      format(stats::median(num(r$s_frm), na.rm = TRUE), digits = 4), "\n")
  tab <- do.call(rbind, lapply(seq_along(truth), function(i) {
    nm <- names(truth)[i]
    e <- num(r[[paste0("err_", nm)]])
    s <- num(r[[paste0("se_", nm)]])
    cov <- abs(e) < stats::qnorm(0.975) * s
    data.frame(coef = nm, term = term[((i - 1L) %% 3L) + 1L],
               truth = truth[[nm]], bias = mean(e),
               empSD = stats::sd(e), meanSE = mean(s),
               se_over_sd = mean(s) / stats::sd(e),
               cover95 = 100 * mean(cov), nNA = sum(is.na(e) | is.na(s)))
  }))
  print(tab, row.names = FALSE, digits = 4)
  k <- sum(num(r$coef_cover_9))
  cat("overall from coef_cover_9 :", mci(k, 9 * nrow(r)), "\n")
  k2 <- sum(round(tab$cover95 / 100 * nrow(r)))
  cat("overall from err_/se_     :", mci(k2, 9 * nrow(r)), "\n")
  cat("the two agree:", k == k2, "\n")
  invisible(tab)
}

ti <- per_block("latent-2p4-recheck.tsv", "TUNING block (in sample)")
to <- per_block("latent-2p4-recheck-oos.tsv", "FRESH block (out of sample)")

hr("Q1 and Q2: the same two coefficients, on both blocks")
cmp <- data.frame(coef = ti$coef, term = ti$term, truth = ti$truth,
                  tune_cov = ti$cover95, tune_se_sd = ti$se_over_sd,
                  fresh_cov = to$cover95, fresh_se_sd = to$se_over_sd)
cmp <- cmp[order(cmp$fresh_cov), ]
print(cmp, row.names = FALSE, digits = 4)
cat("\n`se_over_sd` below 1 means the reported standard error is UNDER\n")
cat("the spread the estimator actually has. Two rows the lane names:\n")
low <- cmp[cmp$fresh_cov < 92, , drop = FALSE]
print(low, row.names = FALSE, digits = 4)
cat("\nthe same two on the TUNING block: coverage",
    paste(format(low$tune_cov, digits = 4), collapse = " and "),
    " se/sd", paste(format(low$tune_se_sd, digits = 4),
                    collapse = " and "), "\n")
cat("\nMonte Carlo half-width on ONE coefficient at 200 replicates,",
    "at p = 0.95:",
    sprintf("%.1f points", 100 * 1.96 * sqrt(0.95 * 0.05 / 200)), "\n")
cat("so a single row at", format(min(low$fresh_cov), digits = 3),
    "percent is",
    format((95 - min(low$fresh_cov)) / (100 * sqrt(0.95 * 0.05 / 200)),
           digits = 3), "standard errors below nominal.\n")

hr("pooling the two blocks, one coefficient at a time")
pool <- do.call(rbind, lapply(seq_along(truth), function(i) {
  nm <- names(truth)[i]
  ri <- rd("latent-2p4-recheck.tsv"); ro <- rd("latent-2p4-recheck-oos.tsv")
  e <- c(num(ri[[paste0("err_", nm)]]), num(ro[[paste0("err_", nm)]]))
  s <- c(num(ri[[paste0("se_", nm)]]), num(ro[[paste0("se_", nm)]]))
  k <- sum(abs(e) < stats::qnorm(0.975) * s)
  p <- k / length(e)
  data.frame(coef = nm, term = term[((i - 1L) %% 3L) + 1L],
             n = length(e), covered = k, cover95 = 100 * p,
             lo = 100 * (p - 1.96 * sqrt(p * (1 - p) / length(e))),
             hi = 100 * (p + 1.96 * sqrt(p * (1 - p) / length(e))),
             se_over_sd = mean(s) / stats::sd(e))
}))
print(pool, row.names = FALSE, digits = 4)
cat("\nrows whose 400-replicate interval EXCLUDES 95:\n")
print(pool[pool$hi < 95, ], row.names = FALSE, digits = 4)
cat("\noverall over 400 replicates:",
    mci(sum(pool$covered), sum(pool$n)), "\n")
