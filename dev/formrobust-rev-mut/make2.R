# Re-check (punch round 1) mutation jobs: the review's 25 mutants on
# their files, the worker's P01 to P06 on theirs, and new mutants R01 to
# R08 for the round-1 code on every file of the lane.
dir <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/"
H <- sprintf('source("%s/helper.R", local = TRUE)', dir)
R <- list(
  R01_epred_fill_off = 'put("arma_cond_fill_epred", function(fit, rspec, newdata, re_formula, dpar = NULL) NULL)',
  R02_mv_grid_switch_off = 'mut("emm_target", "design && length(targets) > 1L", "FALSE && length(targets) > 1L")',
  R03_cens_y2_length_off = 'mut("assemble_frame", "if (!is.null(nd_) && length(v) != nd_) {", "if (FALSE) {")',
  R04_function_var_off = 'mut("expr_frame_vars", "if (is.function(x)) {", "if (FALSE) {")',
  R05_mv_update_allowed = 'mut("update.frmtmb_fit", "frm_stop(\\"Updating formulas of multivariate models is not yet \\",", "if (FALSE) frm_stop(\\"Updating formulas of multivariate models is not yet \\",")',
  R06_call_stores_object = 'put("bform_call", function(b) NULL)',
  R07_bare_const_to_frame = 'mut("assemble_frame", "if (is.name(a)) {\\n for (v in expr_frame_vars(a, data, resp$formula_env)) {\\n add_part(as.name(v))\\n }", "if (is.name(a)) {\\n add_part(a)")',
  R08_family_call_always = 'put("family_call_of", function(fam) fam)'
)
files <- c("tests/testthat/test-arma-na-newdata.R", "tests/testthat/test-aterm-expr.R",
           "tests/testthat/test-bernoulli-coding.R",
           "tests/testthat/test-brms-api-formrobust.R",
           "tests/testthat/test-offset-grid.R", "tests/testthat/test-update-pool.R",
           "tests/testthat/test-open-issues.R",
           "extensions/frmtmb.sample/tests/testthat/test-formrobust-draws.R")
pkg_of <- function(f) if (grepl("frmtmb.sample", f)) "frmtmb.sample" else "frmtmb"
jobs <- character(0)
for (nm in names(R)) {
  writeLines(c(H, R[[nm]]), file.path(dir, paste0(nm, ".R")))
  for (f in files) jobs <- c(jobs, paste(nm, pkg_of(f), paste0(wt, f)))
}
old <- readLines(paste0(wt, "dev/formrobust-rev-log/jobs-mut.txt"))
pj <- readLines(paste0(wt, "dev/formrobust-log/p1-jobs-mut.txt"))
pj <- pj[grepl("^P0", pj)]
pj <- sub(" (tests|extensions)/", paste0(" ", wt, "\\1/"), pj)
for (p in unique(sub(" .*$", "", pj))) {
  file.copy(paste0(wt, "dev/formrobust-p1-mut/", p, ".R"),
            file.path(dir, paste0(p, ".R")), overwrite = TRUE)
}
all <- c(old, pj, jobs)
writeLines(all, paste0(wt, "dev/formrobust-rev-log/jobs-mut2.txt"))
cat(length(all), "jobs\n")
