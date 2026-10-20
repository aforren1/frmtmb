# Mutation check for tests/testthat/helper-compat-contract.R, the check
# test-compat-extensions.R runs in a child with all seven extensions
# loaded. Each mutant plants what the review (dev/reviews/
# 2026-10-14-testiso.md, finding 1 and R1) found the narrowed
# test-compat.R could no longer see, and the check must name it. The
# control (none) holds frmtmb.spline's frm_curve_contrast row, which
# must pass: it is frm_curve()'s contrast argument.
#
#   Rscript dev/testiso-contract-mutants.R <mutant>
.libPaths(c("C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (p in c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent",
            "frmtmb.learn", "frmtmb.ode", "frmtmb.sample",
            "frmtmb.spline")) loadNamespace(p)
mut <- commandArgs(TRUE)[1]
switch(mut,
  none = NULL,
  # a contributed rule whose side is a misspelt expects = name
  typo = frmtmb_register_compat(rules = function() {
    b <- compat_rule_builder()
    b$r("gaussian", "hmmm", "works", "mutant: misspelt expects name")
    b$rules()
  }, expects = "hmmm"),
  # a contributed method row that names no function
  method = frmtmb_register_compat(features = c(zz_contrib_method = "method")),
  # an argument row whose argument is misspelt
  misarg = frmtmb_register_compat(
    features = c(frm_curve_contrsat = "method")),
  # an argument row naming the dots, which every f(...) would accept
  dots = frmtmb_register_compat(features = c("frm_sample_..." = "method")),
  # a contributed rule shadowed by a later one: it wins no pair
  dead = frmtmb_register_compat(rules = function() {
    b <- compat_rule_builder()
    b$r("royston_parmar", "cens()", "untested", "mutant, shadowed")
    b$r("royston_parmar", "cens()", "works", "mutant, wins")
    b$rules()
  }),
  stop("unknown mutant ", mut))
source("tests/testthat/helper-compat-contract.R")
v <- compat_contract_violations()
ft <- frm_compat_features()
cat("CONTRACT", mut, "| spline",
    as.character(packageVersion("frmtmb.spline")),
    basename(dirname(find.package("frmtmb.spline"))),
    "| frm_curve_contrast row:", "frm_curve_contrast" %in% ft$name,
    "| violations:", length(v), "\n")
for (x in v) cat("  ", x, "\n")
