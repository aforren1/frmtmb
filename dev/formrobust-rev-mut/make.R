# Writes one mutant file per mutant into this directory, and the job
# list dev/formrobust-rev-log/jobs-mut.txt (<mutant> <pkg> <test file>).
dir <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut"
tt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/"
st <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/extensions/frmtmb.sample/tests/testthat/"
H <- sprintf('source("%s/helper.R", local = TRUE)', dir)
M <- list(
  M01_expr_vars_unfiltered = list(
    'put("expr_frame_vars", function(a, data, env) if (is.character(a)) a else all.vars(a))',
    "test-aterm-expr.R"),
  M02_no_recycle = list(
    'mut("assemble_frame", "if (length(v) == 1L) v <- rep(v, n)", "NULL")',
    "test-aterm-expr.R"),
  M03_scale_no_length = list(
    'mut("parse_response", "as.call(list(len_, w))", "1")',
    "test-aterm-expr.R"),
  M04_offset_vars_none = list(
    'put("offset_vars", function(f) character(0))',
    c("test-offset-grid.R")),
  M05_emm_keep_offsets = list(
    'put("emm_drop_offsets", function(object) object)',
    "test-offset-grid.R"),
  M05b_emm_terms_keep_offset = list(
    'mut("emm_terms", "attr(tt, \\"offset\\") <- NULL", "NULL")',
    "test-offset-grid.R"),
  M06_fill_expected_in_predict = list(
    'mut("arma_cond_fill_dpars", "autocor_cond_mu_fill(fit, acn, dpv[[\\"mu\\"]], y, draw)", "autocor_cond_mu_fill(fit, acn, dpv[[\\"mu\\"]], y)")',
    "test-arma-na-newdata.R"),
  M07_err_omits_ma = list(
    'mut("autocor_cond_mu_fill", "err[[t]] <- yt - mu[rows] - sma", "err[[t]] <- yt - mu[rows]")',
    "test-arma-na-newdata.R"),
  M07b_fill_at_unshifted_mean = list(
    'mut("autocor_cond_mu_fill", "if (is.null(draw)) m[na]", "if (is.null(draw)) mu[rows][na]")',
    "test-arma-na-newdata.R"),
  M08_bern_levels_reversed = list(
    'old <- get("bernoulli_levels", asNamespace("frmtmb")); put("bernoulli_levels", function(y) rev(old(y)))',
    "test-bernoulli-coding.R"),
  M09_newdata_not_coded = list(
    'put("response_codes_newdata", function(rspec, y, what) y)',
    "test-bernoulli-coding.R"),
  M10_levels_not_carried = list(
    'mut("assemble_frame", "resp$family[[\\"bin_levels\\"]] <- bin_lv", "NULL")',
    "test-bernoulli-coding.R"),
  M10b_carry_dropped = list(
    'mut("carry_finalized_responses", "spec$responses[[rn_]]$family[[\\"bin_levels\\"]] <- fr_$family[[\\"bin_levels\\"]]", "NULL")',
    "test-bernoulli-coding.R"),
  M11_no_pool = list(
    'put("update_complete_bform", function(object, new) NULL)',
    "test-update-pool.R"),
  M12_pool_keeps_old = list(
    'mut("update_pool_pars", "fromLast = TRUE", "fromLast = FALSE")',
    "test-update-pool.R"),
  M13_no_ones_column = list(
    'put("rsv_lower_fill", function(spec, data) data)',
    "test-brms-api-formrobust.R"),
  M14_no_newdata_ones = list(
    'mut("pred_design", "newdata[[\\"intercept\\"]] <- rep(1, nrow(newdata))", "NULL")',
    "test-brms-api-formrobust.R"),
  M15_autocor_arg_ignored = list(
    'put("add_ac_terms", function(f, autocor) f)',
    "test-brms-api-formrobust.R"),
  M16_autocor_no_warning = list(
    'mut("autocor.frmtmb_fit", "frm_warning(", "invisible(")',
    "test-brms-api-formrobust.R"),
  M17_offset_var_displayed = list(
    'mut("ce_lp_vars", "if (!rsv && length(attr(tt, \\"offset\\")))", "if (FALSE)")',
    "test-offset-grid.R"),
  M18_dul_ignored = list(
    'mut("assemble_frame", "drop.unused.levels = drop_unused_levels", "drop.unused.levels = TRUE", times = 2L)',
    "test-brms-api-formrobust.R"),
  M19_bern_fraction_guard_off = list(
    'mut("extract_y", "any(y != round(y), na.rm = TRUE)", "FALSE")',
    c("test-bernoulli-coding.R", "test-open-issues.R")),
  M20_ac_term_added_allowed = list(
    'put("stop_ac_term_added", function(what = "brmsformula") invisible(NULL))',
    "test-brms-api-formrobust.R"),
  S01_draws_fill_expected = list(
    'put("arma_cond_fill_dpars", function(fit, rspec, newdata, dpars_fn) dpars_fn(fit))',
    "S:test-formrobust-draws.R"),
  S02_draws_newdata_not_coded = list(
    'put("response_codes_newdata", function(rspec, y, what) y)',
    "S:test-formrobust-draws.R")
)
jobs <- character(0)
for (nm in names(M)) {
  f <- file.path(dir, paste0(nm, ".R"))
  writeLines(c(H, M[[nm]][[1]]), f)
  for (tf in M[[nm]][[2]]) {
    if (startsWith(tf, "S:")) {
      jobs <- c(jobs, paste(nm, "frmtmb.sample", paste0(st, sub("^S:", "", tf))))
    } else {
      jobs <- c(jobs, paste(nm, "frmtmb", paste0(tt, tf)))
    }
  }
}
writeLines(jobs, "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-log/jobs-mut.txt")
cat(length(M), "mutants,", length(jobs), "jobs\n")
