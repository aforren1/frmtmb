# Survey of bin 2P: the 270 assertions of brms 2.23.0's suite that read a
# brms object frmtmb builds differently. Groups them by the object they
# read, says what adapter frmtmb would need, and classifies every
# assertion that is not ported.
#
#   Rscript dev/brmsport-bin2-survey.R
#
# The names each block calls are READ from its parsed source, not typed:
# every function name in the block that brms 2.23.0 defines is listed
# as exported or internal. Writes dev/brmsport-bin2-groups.tsv (one row
# per group) and dev/brmsport-bin2-survey.tsv (one row per assertion),
# which dev/brmsport-bin2-ledger.R reads.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
source("dev/brmsport-bin2-blocks.R")
source("dev/brmsport-bin2-adapters-spec.R")
topic_of <- function(f) sub("^tests[.](.*)[.]R$", "\\1", f)
bns <- asNamespace("brms")
bexp <- getNamespaceExports("brms")

blocks <- brms_bin_blocks("2P")
group_of <- function(b) {
  f <- b$file; l <- b$label
  if (f == "tests.distributions.R") {
    if (grepl("ordinal|categorical|sratio|cratio", l)) return("G2")
    return("G1")
  }
  switch(f,
         "tests.log_lik.R" = "G3",
         "tests.posterior_predict.R" = "G4",
         "tests.brmsterms.R" = "G5",
         "tests.brmsfit-helpers.R" =
           if (grepl("autocorrelation", l)) "G6" else "G7",
         stop("no group for ", f, ": ", l))
}

# What each group reads, what frmtmb has that it could map onto, the
# cost, and the verdict class of an assertion left unported.
groups <- data.frame(stringsAsFactors = FALSE, rbind(
  c("G1",
      paste0(
        "brms's exported d/p/q/r functions of its own families ",
        "(dstudent_t, dskew_normal, dvon_mises, dbeta_binomial, ",
        "dzero_inflated_*, dhurdle_*, dwiener, ...) and the ",
        "multivariate normal and Student t"),
      paste0(
        "every exported family constructor's $lpdf(y, dpars, ",
        "aterms) (a log density) and $sim(dpars, aterms, n); ",
        "$lcdf where the family has one. No quantile function, ",
        "and no multivariate density outside a fit"),
      paste0(
        "low for the d and r functions of the families frmtmb has: ",
        "a shim per brms signature over row_lpdf() and sim_draw() ",
        "(ported, fix round 2). p needs a CDF, which only ",
        "inverse.gaussian has; q has no counterpart"),
      # the unported G1 rows are the multivariate normal and t block
      # (this class) and the wiener block (its own override below)
      "no-density",
      paste0(
        "decision: frmtmb has no multivariate normal or t density ",
        "outside a fit (the rescor density is built inside the ",
        "objective), and SPEC.md:153-154 puts numeric d/p/q/r ",
        "functions outside frmtmb, in gamlss.dist or extraDistr; ",
        "brms's 8 validation messages of these functions have no ",
        "counterpart either")),
  c("G2",
      paste0(
        "brms's unexported ordinal and categorical densities and ",
        "link maps on arrays (brms:::dcumulative, dsratio, ",
        "link_cumulative, inv_link_sratio, dcategorical, ",
        "link_categorical, ...), called unqualified inside ",
        "brms's own test environment"),
      paste0(
        "the ordinal densities: each ordinal family object's ",
        "exported ord_lpdf_make(identity), a thin read (ported). ",
        "The link maps on arrays and the categorical functions: ",
        "none, frmtmb exports no map between linear predictors ",
        "and category probabilities"),
      paste0(
        "low for the two density blocks; high for the link maps ",
        "(the array convention, draws x observations x ",
        "categories, would be new frmtmb functions)"),
      "brms-internal",
      paste0(
        "reads brms's unexported ordinal and categorical array ",
        "densities and link maps; frmtmb has no counterpart")),
  c("G3",
      paste0(
        "brms:::log_lik_<family>(i, prep), with prep a ",
        "hand-built brmsprep (dpars as draws x observations ",
        "matrices, data$Y, trials, weights, cens)"),
      paste0(
        "every exported family constructor's $lpdf(y, dpars, ",
        "aterms), documented in ?frmtmb_family: a thin adapter ",
        "evaluates it per draw at observation i; an ordinal ",
        "family's ord_lpdf_make(identity) takes brms's ",
        "thresholds. Autocorrelation (ARMA, SAR, FCOR), ",
        "multivariate, categorical and custom-family preps have ",
        "no such read"),
      paste0(
        "low for the univariate families (one adapter, a dpar ",
        "name map and brms's reference densities bound as data); ",
        "high for the rest"),
      "not-ported",
      paste0(
        "not ported, by cost: the prep needs more than the thin ",
        "adapter reads (ARMA, SAR or FCOR structure, a multivariate ",
        "response, wiener's frmtmb.eam family, a custom family), and ",
        "frmtmb has each of these models, so the row is transferable ",
        "and stays in the rate's denominator")),
  c("G4",
      "brms:::posterior_predict_<family>(i, prep), the same brmsprep",
      paste0(
        "every exported family constructor's $sim(dpars, aterms, ",
        "n), documented in ?frmtmb_family: a thin adapter draws ",
        "per draw at observation i; an ordinal family's ",
        "ord_sim_make(identity). The same exclusions as G3"),
      "low for the univariate families; high for the rest",
      "not-ported",
      paste0(
        "not ported, by cost, as in G3: frmtmb simulates each of ",
        "these models, so the row is transferable and stays in the ",
        "rate's denominator")),
  c("G5",
      paste0(
        "brms:::brmsterms(), check_re_formula(), ",
        "update_re_terms(): brms's parsed-formula object and its ",
        "formula editors"),
      paste0(
        "none exported: frm(dry_run = 'frame') returns frmtmb's ",
        "assembled frame, whose structure is not brmsterms's (no ",
        "allvars, no dpars list in brms's shape); re_form_arg() ",
        "resolves a re_formula against a fit, not against a ",
        "formula"),
      "high: an adapter would re-derive brms's object from frmtmb's frame",
      "brms-internal",
      "reads brms:::brmsterms() or its formula editors, brms internals"),
  c("G6",
      paste0(
        "brms:::get_cor_matrix_ar1(), _ma1(), _arma1(), _cosy(), ",
        "_ident(): autocorrelation matrices from parameter draws"),
      paste0(
        "autocor_matrix(fit) is exported but reads a FITTED ",
        "model's correlation, and returns a correlation matrix ",
        "where brms's AR and ARMA matrices are the covariance at ",
        "unit innovation variance; there is no function of the ",
        "parameter"),
      paste0(
        "medium: an adapter would fit a model at the given ",
        "parameter or call internals"),
      "brms-internal",
      paste0(
        "reads brms:::get_cor_matrix_*(), brms internals; ",
        "autocor_matrix() takes a fit")),
  c("G7",
      paste0(
        "brms:::insert_refcat(): the reference category of a ",
        "categorical linear predictor"),
      "none",
      "high: no counterpart",
      "brms-internal",
      "reads brms:::insert_refcat(), a brms internal")))
names(groups) <- c("group", "reads", "frmtmb_counterpart", "cost",
                   "class", "reason")

rows <- list()
for (b in blocks) {
  g <- group_of(b)
  nm <- unique(unlist(lapply(b$stmts, all.names)))
  bn <- nm[vapply(nm, function(n) {
    exists(n, envir = bns, inherits = FALSE, mode = "function")
  }, NA)]
  bn <- setdiff(bn, c("expect_equal", "expect_true", "expect_error",
                      "expect_length", "expect_warning", "expect_is",
                      "SW", "SM"))
  ported <- b$file %in% names(bin2_adapted) &&
    b$label %in% bin2_adapted[[b$file]]$blocks
  for (j in seq_along(b$stmts)) {
    # a statement can hold several assertions (a loop, SW(expect_...));
    # it is one ledger row carrying their count, as the audit counts them
    k <- n_expect(b$stmts[[j]])
    if (k == 0L) next
    # A block whose unported rows take a class other than their group's.
    # brms's G1 wiener block compares brms's two wiener backends
    # (RWiener, rtdists); frmtmb.eam has one implementation, so a shim
    # would make the two calls the same one and the assertion hollow.
    # frmtmb.eam's density agrees with both backends to 5.7e-14 on the
    # log scale at the block's inputs (dev/remeasure-bin2-wiener-g1.R).
    ov <- if (grepl("^wiener distribution", b$label)) {
      c("brms-backends", paste0(
        "decision: the block compares brms's two wiener backends with ",
        "each other; frmtmb.eam has one implementation, which agrees ",
        "with both to 5.7e-14 on the log scale at the block's inputs ",
        "(dev/remeasure-bin2-wiener-g1.R)"))
    } else if (b$label %in% c("dcategorical() works correctly",
               "inv_link_categorical() works correctly for arrays")) {
      # fix 2b: the categorical probability map, which frmtmb's
      # categorical family has; the same class as the categorical
      # log_lik and posterior_predict rows
      c("not-ported", paste0(
        "not ported (cost): brms's categorical probability map, which ",
        "frmtmb's categorical family computes; the existing adapter ",
        "refuses brms's mu1/mu2 layout, as for the categorical log_lik ",
        "and posterior_predict rows"))
    } else {
      c("", "")
    }
    rows[[length(rows) + 1L]] <- data.frame(
      id = sprintf("%s:%d", topic_of(b$file), b$lines[j]), file = b$file,
      n = k, block = gsub("[\r\n\t]+", " ", b$label), group = g,
      ported = ported, class = ov[1], reason = ov[2],
      brms_exported = paste(intersect(bn, bexp), collapse = " "),
      brms_internal = paste(setdiff(bn, bexp), collapse = " "),
      stringsAsFactors = FALSE)
  }
}
S <- do.call(rbind, rows)
stopifnot(sum(S$n) == 270L, !anyDuplicated(S$id))
gs <- tapply(S$n, S$group, sum)
groups$assertions <- as.integer(gs[groups$group])
groups$blocks <- vapply(groups$group, function(g) {
  length(unique(S$block[S$group == g]))
}, 1L)
groups$ported <- vapply(groups$group, function(g) {
  sum(S$n[S$ported & S$group == g])
}, 1L)
utils::write.table(S, "dev/brmsport-bin2-survey.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
utils::write.table(groups, "dev/brmsport-bin2-groups.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("2P assertions:", sum(S$n), "in", nrow(S), "statements of",
    length(blocks), "blocks; ported:", sum(S$n[S$ported]), "\n\n")
print(groups[c("group", "blocks", "assertions", "ported", "class")],
      row.names = FALSE)
cat("\nbrms functions read, per group (from the parsed blocks):\n")
for (g in groups$group) {
  s <- S[S$group == g, ]
  cat(sprintf("%s exported: %s\n   internal: %s\n", g,
              paste(sort(unique(unlist(strsplit(s$brms_exported, " ")))),
                    collapse = " "),
              paste(sort(unique(unlist(strsplit(s$brms_internal, " ")))),
                    collapse = " ")))
}
