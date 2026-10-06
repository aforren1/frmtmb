# Counts for the findings from dev/nanse-mo-lane.R's TSVs: the base
# build and the lane build, seeds 1 to 200.
#   Rscript dev/nanse-mo-lane-sum.R <base.tsv> <lane.tsv>
a <- commandArgs(TRUE)
B <- utils::read.delim(a[1], stringsAsFactors = FALSE)
L <- utils::read.delim(a[2], stringsAsFactors = FALSE)
for (fm in c("int", "main")) {
  b <- B[B$form == fm, ]
  l <- L[L$form == fm, ]
  stopifnot(identical(b$seed, l$seed))
  cat(sprintf("\n== %s: %d seeds\n", fm, nrow(l)))
  for (nm in c("base", "lane")) {
    x <- if (nm == "base") b else l
    alln <- x$se_finite == 0
    some <- x$se_finite > 0 & x$se_finite < x$n_par
    cat(sprintf(paste0("%s: all SEs non-finite %d, some %d, all finite %d;",
                       " SE warning at fit %d; any fit warning %d;",
                       " conditional_effects print fails %d\n"),
                nm, sum(alln), sum(some), sum(x$se_finite == x$n_par),
                sum(x$se_warn), sum(x$n_warn_fit > 0), sum(!x$ce_ok)))
  }
  lost <- l$se_finite < l$n_par
  cat(sprintf("lane: SE warning on seeds with a lost SE: %d of %d;",
              sum(l$se_warn & lost), sum(lost)))
  cat(sprintf(" on seeds with every SE finite: %d of %d\n",
              sum(l$se_warn & !lost), sum(!lost)))
  cat(sprintf("lane: lost seeds whose fit already warned otherwise: %d\n",
              sum(lost & !l$se_warn & l$n_warn_fit > 0)))
  cat(sprintf("lane: lost seeds with no warning at all: %d\n",
              sum(lost & l$n_warn_fit == 0)))
  if (any(lost)) {
    cat("lane lost sets:\n")
    print(table(l$lost[lost]))
  }
  cat("lane: other fit warnings on lost seeds:\n")
  ow <- l$warn_fit[lost & !l$se_warn & l$n_warn_fit > 0]
  print(table(substr(ow, 1, 60)))
  fin <- l$se_finite == l$n_par
  cat(sprintf(paste0("lane: max |SE/ref - 1| over the coefficients and ",
                     "sigma, seeds with every SE finite: %.3g (median %.3g)",
                     "\n"), max(l$max_rel_b[fin], na.rm = TRUE),
              stats::median(l$max_rel_b[fin], na.rm = TRUE)))
  cat("lane: seeds whose base SEs were all finite and now differ:",
      sum(b$se_finite == b$n_par & b$se_b != l$se_b), "\n")
  cat("base->lane change of the all-NaN seeds:\n")
  print(table(base = ifelse(b$se_finite == 0, "allNaN",
                            ifelse(b$se_finite < b$n_par, "some", "finite")),
              lane = ifelse(l$se_finite == 0, "allNaN",
                            ifelse(l$se_finite < l$n_par, "some", "finite"))))
}
