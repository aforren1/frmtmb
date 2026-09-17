## Punch round 1, MINOR 1: the pins must fail when the names lie. Apply
## one mutation in memory (assignInNamespace, nothing installed), then
## run the pinning test files and report PASS and FAIL per file.
##   Rscript dev/brmsnames-mutants.R <mutant> <testfile>...
## Mutants (after dev/brmsnames-rev-mutants.R, adapted to the current
## code, plus the guards the review found unpinned):
##   none      no mutation
##   rlev      brms_r_labels() writes the levels in reverse order
##   rcoef     brms_r_labels() writes the coefficient names in reverse
##   bcoef     brms_coef_table() names a predictor's b_ columns in reverse
##   corord    draws_cor_index() returns R's column-major lower triangle
##   natsig    draws_to_natural() leaves sigma on the log scale
##   residse   VarCorr() on a fit drops betad from the delta method, so
##             residual__ has Est.Error 0
##   llpw      log_lik() stops refusing pointwise = TRUE
##   norename  hypothesis() parses without brms's renaming
## Added when round 1 was verified, for the items that had no mutant:
##   charmap   brms_rename() keeps + - * / ^ = $ (MAJOR 1)
##   stanname  brms_stan_name() keeps a response's _ and . (MAJOR 1)
##   lvldot    r_ levels keep their whitespace (MAJOR 1)
##   lvljoin   an interaction group's levels keep their : (MAJOR 1)
##   nodupref  brms_rename(check_dup = TRUE) stops refusing (MAJOR 2)
##   nosuffix  no __1 suffix on a name two predictors share (MAJOR 2)
##   noredup   brms_check_re_dups() stops refusing (MAJOR 2)
##   natname   an unmodeled sigma is named b_sigma_Intercept (MAJOR 3)
##   rrfill    ranef() leaves a reduced-rank block unfilled (MAJOR 4)
##   lapref    ranef() stops refusing laplace draws by name (MAJOR 4)
##   mvresid   no residual__ on a multivariate VarCorr() (MINOR 4)
##   mvr2      bayes_R2() answers the first response only (MINOR 4)
##   barenp    frm_simulate(newparams =) maps a bare name (USER)
## Added in punch round 2, each restoring the round-1 behavior of one fix:
##   mixnat    a mixture's theta<k> is elementwise natural again
##   hmmnat    an hmm() transition logit is flagged natural again
##   rdup      draws labels lose brms's __1 suffix
##   lvlmerge  merged interaction-group levels are not refused
##   respdup   responses brms spells alike are not refused
##   hypdup    draws hypothesis() reads the first of two same-named columns
##   sdscolon  a by-smooth's sds_ name keeps its ':'
##   bspname   a monotonic scale is b_ again
av <- commandArgs(trailingOnly = TRUE)
mut <- av[1]
files <- av[-1]
source("dev/brmsnames-libs.R")
brmsnames_libs("lane")
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(testthat))
q(library(frmtmb))
q(library(frmtmb.sample))

# Replace one line of a function's deparsed body, and refuse a pattern
# that does not match exactly once, so a mutant cannot silently be a
# no-op
mutate <- function(fun, pkg, from, to) {
  ns <- asNamespace(pkg)
  src <- deparse(get(fun, ns), width.cutoff = 500L)
  hit <- grep(from, src, fixed = TRUE)
  if (length(hit) != 1L) {
    stop("mutant ", mut, ": pattern matched ", length(hit), " lines")
  }
  src[hit] <- sub(from, to, src[hit], fixed = TRUE)
  g <- eval(parse(text = src))
  environment(g) <- ns
  assignInNamespace(fun, g, pkg)
  # frmtmb.sample holds its own copy of an imported function
  if (pkg == "frmtmb" && exists(fun, parent.env(asNamespace("frmtmb.sample")),
                                inherits = FALSE)) {
    imp <- parent.env(asNamespace("frmtmb.sample"))
    unlockBinding(fun, imp)
    assign(fun, g, imp)
  }
}

invisible(switch(mut,
  none = NULL,
  rlev = mutate("brms_r_labels", "frmtmb",
                "lev <- brms_levels(bk, for_r = TRUE)",
                "lev <- rev(brms_levels(bk, for_r = TRUE))"),
  rcoef = mutate("brms_r_labels", "frmtmb",
                 "parts$coef[ci]", "rev(parts$coef)[ci]"),
  bcoef = mutate("brms_coef_table", "frmtmb",
                 "out <- paste0(\"b_\", brms_usc(pre, brms_rename(cn)))",
                 "out <- paste0(\"b_\", brms_usc(pre, brms_rename(rev(cn))))"),
  corord = mutate("draws_cor_index", "frmtmb.sample",
                  "for (i in seq_len(K)[-1]) {",
                  "return(which(lower.tri(diag(K)))); for (i in seq_len(K)[-1]) {"),
  natsig = mutate("draws_to_natural", "frmtmb.sample",
                  "nc <- draws_natural_cols(fit)", "return(m)"),
  residse = mutate("VarCorr.frmtmb_fit", "frmtmb",
                   "c(\"theta\", \"betad\", \"thetar\")",
                   "c(\"theta\", \"thetar\")"),
  llpw = mutate("log_lik.frmtmb_draws", "frmtmb.sample",
                "if (pointwise) {", "if (FALSE) {"),
  norename = mutate("hyp_rename", "frmtmb",
                    "brms_rename(x, c(\":\", \"[\", \"]\", \",\"), c(\"___\", \".\", \".\", \"..\"))",
                    "x"),
  charmap = mutate("brms_rename", "frmtmb",
                   "c(rep(\"\", 9), \"P\", \"M\", \"MU\", \"D\", \"E\", \"EQ\", \"USD\")",
                   "c(rep(\"\", 9), \"+\", \"-\", \"*\", \"/\", \"^\", \"=\", \"$\")"),
  stanname = mutate("brms_stan_name", "frmtmb",
                    "gsub(\"[._]\", \"\", make.names(x, unique = TRUE))",
                    "make.names(x, unique = TRUE)"),
  lvldot = mutate("brms_levels", "frmtmb", "if (for_r) ", "if (FALSE) "),
  lvljoin = mutate("brms_levels", "frmtmb",
                   "lev <- gsub(\":\", \"_\", lev, fixed = TRUE)",
                   "lev <- lev"),
  nodupref = mutate("brms_rename", "frmtmb",
                    "if (check_dup && any(dup)) {", "if (FALSE) {"),
  nosuffix = mutate("brms_coef_table", "frmtmb",
                    "make.unique(c(nm$beta, nm$betad[keep_d]), sep = \"__\")",
                    "c(nm$beta, nm$betad[keep_d])"),
  noredup = mutate("brms_check_re_dups", "frmtmb",
                   "seen <- character(0)", "return(invisible(NULL))"),
  natname = mutate("brms_dpar_written", "frmtmb",
                   "if (identical(lp[[\"par\"]], \"beta\")) ",
                   "if (TRUE) "),
  rrfill = mutate("draws_ranef_fill", "frmtmb.sample",
                  "if (!length(miss)) ", "if (TRUE) "),
  lapref = mutate("draws_ranef_fill", "frmtmb.sample",
                  "if (draws_is_laplace(x)) {", "if (FALSE) {"),
  mvresid = mutate("varcorr_residual_layout", "frmtmb",
                   "if (!mv) {", "if (mv) return(NULL); if (!mv) {"),
  # 0.6.0's shape: one R2 column, whichever response came first
  mvr2 = mutate("bayes_R2.frmtmb_draws", "frmtmb.sample",
                "R2 <- lapply(sel, function(r) {",
                "sel <- sel[1L]; R2 <- lapply(sel, function(r) {"),
  # the bare spelling mapped to brms's before anything reads the names,
  # which is what accepting both spellings would take
  barenp = mutate("frm_simulate", "frmtmb",
                  "check_natural_supported(frame)",
                  paste0("check_natural_supported(frame); ",
                         "lg <- attr(slots, \"legacy\"); ",
                         "hit <- names(np_natural) %in% names(lg); ",
                         "names(np_natural)[hit] <- ",
                         "lg[names(np_natural)[hit]]")),
  mixnat = mutate("brms_coef_table", "frmtmb",
                  "if (lp[[\"dpar\"]] %in% mix && grepl(",
                  "if (FALSE && grepl("),
  hmmnat = mutate("brms_coef_table", "frmtmb",
                  "!lp[[\"dpar\"]] %in% fam[[\"link_scale_dpars\"]]",
                  "TRUE"),
  rdup = mutate("brms_par_labels", "frmtmb",
                "make.unique(out, sep = \"__\")", "out"),
  lvlmerge = mutate("brms_check_re_dups", "frmtmb",
                    "if (anyDuplicated(lev)) {", "if (FALSE) {"),
  respdup = mutate("brms_check_re_dups", "frmtmb",
                   "if (length(rs) > 1L && anyDuplicated(brms_stan_name(rs))) {",
                   "if (FALSE) {"),
  hypdup = mutate("hypothesis.frmtmb_draws", "frmtmb.sample",
                  "if (length(twice)) {", "if (FALSE) {"),
  sdscolon = mutate("brms_sds_names", "frmtmb",
                    "gsub(\":\", \"\", brms_rename(",
                    "(function(a, b, x, fixed) x)(\":\", \"\", brms_rename("),
  bspname = mutate("brms_coef_table", "frmtmb",
                   "out[mo] <- paste0(\"bsp_\"", "out[mo] <- paste0(\"b_\""),
  stop("unknown mutant ", mut)))
cat("mutant:", mut, "\n")

for (file in files) {
  pkg <- regmatches(file, regexpr("frmtmb[.][a-z]+", file))
  if (!length(pkg)) pkg <- "frmtmb"
  q(library(pkg, character.only = TRUE))
  res <- q(test_file(file, reporter = "silent", package = pkg))
  df <- as.data.frame(res)
  cat(sprintf("%-22s PASS %4d FAIL %3d ERROR %2d\n", basename(file),
              sum(df$passed), sum(df$failed), sum(df$error)))
  # the message says whether a mutant failed on behavior or only on a
  # crash; an ERROR row alone cannot tell the two apart
  for (i in which(df$failed > 0 | df$error)) {
    cat("    failing:", df$test[i], "\n")
    for (r in res[[i]]$results) {
      if (inherits(r, "expectation_failure") ||
          inherits(r, "expectation_error")) {
        cat("      msg:", substr(gsub("[\r\n]+", " ", conditionMessage(r)),
                                 1, 200), "\n")
      }
    }
  }
}
cat("DONE\n")
