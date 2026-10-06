# Summarize dev/nanse-mo-sweep.R's TSV into the counts the findings
# quote.  Rscript dev/nanse-mo-sweep-sum.R <tsv>
f <- commandArgs(TRUE)[1]
X <- utils::read.delim(f, stringsAsFactors = FALSE)
for (fm in unique(X$form)) {
  Y <- X[X$form == fm, ]
  all_nan <- Y$se_plain_finite == 0
  cat(sprintf("\n== form %s: %d seeds\n", fm, nrow(Y)))
  cat(sprintf("all SEs non-finite (sdreport): %d; some: %d\n",
              sum(all_nan), sum(Y$se_plain_finite > 0 &
                                  Y$se_plain_finite < Y$n_par)))
  cat(sprintf("  of the all-NaN: code 0 %d, no warning at fit %d, both %d\n",
              sum(all_nan & Y$code == 0), sum(all_nan & Y$n_warn == 0),
              sum(all_nan & Y$code == 0 & Y$n_warn == 0)))
  cat(sprintf("  pdHess TRUE on all-NaN: %d\n", sum(all_nan & Y$pdHess)))
  cat(sprintf("seeds with a simplex component < 1e-8: %d; of all-NaN: %d\n",
              sum(Y$n_sat > 0), sum(all_nan & Y$n_sat > 0)))
  cat(sprintf("finite-SE seeds with a saturated component: %d\n",
              sum(!all_nan & Y$n_sat > 0)))
  cat(sprintf("scaled inverse: all finite on %d of the all-NaN seeds\n",
              sum(all_nan & Y$se_scaled_finite == Y$n_par)))
  cat(sprintf("held inverse: SE finite off the held set on %d of %d all-NaN\n",
              sum(all_nan & Y$se_held_finite == Y$n_par - Y$n_held),
              sum(all_nan)))
  cat(sprintf("held set non-empty: %d seeds (all-NaN %d, finite %d)\n",
              sum(Y$n_held > 0), sum(all_nan & Y$n_held > 0),
              sum(!all_nan & Y$n_held > 0)))
  cat("held sets on all-NaN seeds:\n")
  print(table(Y$held[all_nan], useNA = "ifany"))
  cat("held sets on finite seeds:\n")
  print(table(Y$held[!all_nan], useNA = "ifany"))
  cat("zero diagonals:", sum(Y$n_zero_diag > 0), "seeds\n")
  cat("rcond(H) on all-NaN seeds: max", format(max(Y$rcond_H[all_nan])),
      "; on finite seeds: min", format(min(Y$rcond_H[!all_nan])), "\n")
  cat("max |gradient| on all-NaN seeds:",
      format(max(Y$maxgrad[all_nan])), "\n")
  cat("all-NaN seeds:", paste(Y$seed[all_nan], collapse = " "), "\n")
}
