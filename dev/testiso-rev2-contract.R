# Reviewer scratch (fix round 1): the reviewer's mutants against
# compat_contract_violations() (tests/testthat/helper-compat-contract.R),
# with all seven extensions loaded, as the child of
# test-compat-extensions.R runs it.
#   Rscript dev/testiso-rev2-contract.R <mutant>
.libPaths(c("C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mut <- commandArgs(TRUE)[1]
suppressMessages(library(frmtmb))
for (p in c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
            "frmtmb.ode", "frmtmb.sample", "frmtmb.spline")) loadNamespace(p)
source("C:/Users/adf44/source/r/frmtmb-wt-testiso/tests/testthat/helper-compat-contract.R")
meth <- function(nm) frmtmb_register_compat(features = setNames("method", nm))
switch(mut,
  none = NULL,
  typo = frmtmb_register_compat(rules = function() {
    b <- compat_rule_builder()
    b$r("gaussian", "hmmm", "works", "mutant: misspelt expects name")
    b$rules()
  }, expects = "hmmm"),
  method = meth("zz_contrib_method"),
  # a misspelt argument of a real function
  misarg = meth("frm_curve(contrsat =)"),
  # the dots of a function, if it has them, are not an argument setting
  dots = meth("frm_curve(... =)"),
  # an exported NON-function of a loaded extension, if any, must not pass
  # as a method; a core function named by an extension row is the other
  # direction (frm_lp_basis is core's)
  corefn = meth("frm_lp_basis2"),
  # a prefix that is an exported function and a suffix that happens to
  # be one of its formals, though nobody meant an argument setting
  coincide = meth("frm_sample_control"),
  coincide2 = meth("frm_curve_seed"),
  # a core function as the prefix: core is not searched
  corepfx = meth("frm_compat_status"),
  stop("unknown mutant ", mut))
cat("frm_curve formals:", paste(names(formals(frmtmb.spline::frm_curve)),
                                collapse = " "), "\n")
v <- compat_contract_violations()
cat("CONTRACT", mut, "| violations:", length(v), "\n")
for (x in v) cat("  ", x, "\n")
