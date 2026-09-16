# Assigns a bin to every test_that block in brms 2.23.0's own suite and
# emits the counts the audit quotes. The counts are GENERATED here and
# pasted verbatim, never typed, because a hand-typed total in a sentence
# is what no verifier checks.
#
# Bins, as Phase 2.6a defines them:
#   1  transfers as-is: formula grammar, priors, family and data
#      validation, argument names and return shapes
#   2  transfers with an adaptation. Two reasons, tracked apart:
#      2T the assertion reads a POSTERIOR that ML replaces with an
#         estimate and a standard error, so it is restated
#      2P the assertion reaches a brms INTERNAL object (brmsprep,
#         brmsterms) that frmtmb implements differently, so it is
#         rewired; its tolerance is unchanged, usually an identity
#   3  does not transfer. Three reasons:
#      3S Stan code or Stan data
#      3M MCMC diagnostics, or a quantity defined only over draws
#      3I a brms-internal utility with no frmtmb counterpart
#
# Tier says which package the block can run against: core frmtmb, the
# frmtmb.sample draws tier, or both.
dir <- "dev/brms-suite/brms/tests/testthat"
files <- sort(list.files(dir, pattern = "[.]R$", full.names = TRUE))

blocks <- list()
for (f in files) {
  ex <- parse(f, keep.source = FALSE)
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    if (!is.call(e) || !identical(e[[1]], as.name("test_that"))) next
    nms <- all.names(e, functions = TRUE, unique = FALSE)
    v <- grep("^expect_[a-z0-9_]+$", nms, value = TRUE)
    v <- v[!grepl("^expected", v)]
    blocks[[length(blocks) + 1L]] <- data.frame(
      file = basename(f),
      desc = gsub("[\r\n\t]+", " ", paste(as.character(e[[2]]),
                                          collapse = " ")),
      expects = length(v),
      refusals = sum(v %in% c("expect_error", "expect_warning",
                              "expect_message")),
      stringsAsFactors = FALSE
    )
  }
}
b <- do.call(rbind, blocks)

# file-level rules; a per-block rule below overrides them
file_rule <- c(
  "tests.brmsformula.R"        = "1|both",
  "tests.families.R"           = "1|both",
  "tests.data-helpers.R"       = "1|core",
  "tests.emmeans.R"            = "1|core",
  "tests.misc.R"               = "3I|-",
  "tests.stop2.R"              = "3I|-",
  "tests.rename_pars.R"        = "3I|-",
  "tests.restructure.R"        = "3I|-",
  "tests.exclude_pars.R"       = "3S|-",
  "tests.stan_functions.R"     = "3S|-",
  "tests.read_csv_as_stanfit.R" = "3S|-",
  "tests.priorsense.R"         = "3M|-",
  "tests.log_lik.R"            = "2P|core",
  "tests.posterior_predict.R"  = "2P|both",
  "tests.distributions.R"      = "2P|core",
  "tests.posterior_epred.R"    = "3I|-",
  "tests.stancode.R"           = "3S|-"
)

# block-level rules: file, regex on the description, bin|tier
blk_rule <- list(
  c("tests.brm.R", "mock backend", "3S|-"),
  c("tests.brm.R", "expected errors", "1|core"),
  c("tests.brmsterms.R", "fixed auxiliary parameters", "1|both"),
  c("tests.brmsterms.R", "unused variables", "1|core"),
  c("tests.brmsterms.R", ".", "2P|core"),
  c("tests.priors.R", "overall intercept priors", "2T|core"),
  c("tests.priors.R", ".", "1|core"),
  c("tests.brmsfit-helpers.R", "probit", "1|core"),
  c("tests.brmsfit-helpers.R", "make_conditions", "1|core"),
  c("tests.brmsfit-helpers.R", "autocorrelation matrices", "2P|core"),
  c("tests.brmsfit-helpers.R", "insert_refcat", "2P|core"),
  c("tests.brmsfit-helpers.R", "evidence_ratio", "3M|-"),
  c("tests.brmsfit-helpers.R", ".", "3I|-"),
  # tests.standata.R
  c("tests.standata.R", "rejects incorrect response", "1|core"),
  c("tests.standata.R", "accepts correct response", "1|core"),
  c("tests.standata.R", "suggests using family bernoulli", "1|core"),
  c("tests.standata.R", "rejects incorrect addition", "1|core"),
  c("tests.standata.R", "correct values for addition terms", "1|core"),
  c("tests.standata.R", "removes NAs", "1|core"),
  c("tests.standata.R", "initial data order", "1|core"),
  c("tests.standata.R", "'poly' function", "1|core"),
  c("tests.standata.R", "Cell-mean coding", "1|core"),
  c("tests.standata.R", "dots in formula", "1|core"),
  c("tests.standata.R", "reserved variables", "1|core"),
  c("tests.standata.R", "unused interval censoring", "1|core"),
  c("tests.standata.R", "drop_unused_factor", "1|core"),
  c("tests.standata.R", "fixed distributional parameters", "1|core"),
  c("tests.standata.R", "'subset' addition argument", "1|core"),
  c("tests.standata.R", "'mi' terms with 'subset'", "1|core"),
  c("tests.standata.R", "grouped ordinal thresholds", "1|core"),
  c("tests.standata.R", "addition term 'rate'", "1|core"),
  c("tests.standata.R", "by variables in grouping terms", "1|core"),
  c("tests.standata.R", ".", "3S|-"),
  # tests.brmsfit-methods.R
  c("tests.brmsfit-methods.R", "^as_draws|^as\\.data\\.frame|^as\\.matrix|^as\\.array|^as\\.mcmc",
    "1|sample"),
  c("tests.brmsfit-methods.R", "^autocor|^conditional_effects has|^conditional_smooths|^family has|^fitted|^fixef|^formula has|^hypothesis|^model\\.frame|^ngrps|^nobs|^predict has|^ranef|^residuals|^variables|^vcov|^contrasts of grouping",
    "1|core"),
  c("tests.brmsfit-methods.R", "^pp_check|^plot has|^prior_summary|^print has|^summary has|^update has",
    "1|core"),
  c("tests.brmsfit-methods.R", "^posterior_summary|^posterior_interval|^posterior_predict|^posterior_linpred|^posterior_epred|^predictive_error|^pp_mixture|^prior_draws|^posterior_samples|^posterior_average|^nsamples|^ndraws and|^pairs has",
    "1|sample"),
  c("tests.brmsfit-methods.R", "^bayes_R2|^coef has|^log_lik has|^VarCorr|^waic|^loo", "2T|both"),
  c("tests.brmsfit-methods.R", "^stancode|^standata|^inits", "3S|-"),
  c("tests.brmsfit-methods.R", "^mcmc_plot|^diagnostic convenience|^bayes_factor|^bridge_sampler|^launch_shinystan|^combine_models|^model_weights|^pp_average|^post_prob",
    "3M|sample"),
  c("tests.brmsfit-methods.R", "^plot of conditional_effects", "3I|-")
)

b$bin <- NA_character_
b$tier <- NA_character_
for (r in blk_rule) {
  hit <- b$file == r[1] & grepl(r[2], b$desc) & is.na(b$bin)
  parts <- strsplit(r[3], "|", fixed = TRUE)[[1]]
  b$bin[hit] <- parts[1]
  b$tier[hit] <- parts[2]
}
for (nm in names(file_rule)) {
  hit <- b$file == nm & is.na(b$bin)
  parts <- strsplit(file_rule[[nm]], "|", fixed = TRUE)[[1]]
  b$bin[hit] <- parts[1]
  b$tier[hit] <- parts[2]
}

# Fail closed. An unclassified block must stop the run rather than be
# dropped from a total, which is how a count silently shrinks.
if (any(is.na(b$bin))) {
  print(b[is.na(b$bin), c("file", "desc")])
  stop("unclassified blocks: ", sum(is.na(b$bin)))
}
write.table(b, "dev/brmssuite-classified.tsv", sep = "\t",
            row.names = FALSE, quote = TRUE)

b$bin1 <- substr(b$bin, 1, 1)
cat("<!-- generated by dev/brmssuite-classify.R; do not edit by hand -->\n")
cat("\n### Per file\n\n")
cat("| file | blocks | assertions | bin 1 | bin 2 | bin 3 | verdict |\n")
cat("|---|---|---|---|---|---|---|\n")
agg <- split(b, b$file)
ord <- names(sort(sapply(agg, function(x) -sum(x$expects))))
for (nm in ord) {
  x <- agg[[nm]]
  v <- paste(sort(unique(x$bin)), collapse = ", ")
  cat(sprintf("| `%s` | %d | %d | %d | %d | %d | %s |\n", nm, nrow(x),
              sum(x$expects), sum(x$expects[x$bin1 == "1"]),
              sum(x$expects[x$bin1 == "2"]),
              sum(x$expects[x$bin1 == "3"]), v))
}
cat(sprintf("| **total** | **%d** | **%d** | **%d** | **%d** | **%d** | |\n",
            nrow(b), sum(b$expects), sum(b$expects[b$bin1 == "1"]),
            sum(b$expects[b$bin1 == "2"]), sum(b$expects[b$bin1 == "3"])))

cat("\n### Per bin\n\n")
cat("| bin | blocks | assertions | share of assertions |\n|---|---|---|---|\n")
for (k in c("1", "2T", "2P", "3S", "3M", "3I")) {
  x <- b[b$bin == k, ]
  cat(sprintf("| %s | %d | %d | %.1f%% |\n", k, nrow(x), sum(x$expects),
              100 * sum(x$expects) / sum(b$expects)))
}
cat("\n### Bin 1 by tier\n\n")
x <- b[b$bin1 == "1", ]
cat("| tier | blocks | assertions |\n|---|---|---|\n")
for (t in sort(unique(x$tier))) {
  y <- x[x$tier == t, ]
  cat(sprintf("| %s | %d | %d |\n", t, nrow(y), sum(y$expects)))
}
cat("\n### Refusal assertions (expect_error/warning/message) by bin\n\n")
cat("| bin | refusal assertions |\n|---|---|\n")
for (k in c("1", "2T", "2P", "3S", "3M", "3I")) {
  cat(sprintf("| %s | %d |\n", k, sum(b$refusals[b$bin == k])))
}
