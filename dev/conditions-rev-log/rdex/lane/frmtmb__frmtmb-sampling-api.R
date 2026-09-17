.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-sampling-api
### Title: Core internals for a sampling extension
### Aliases: frmtmb-sampling-api build_objective row_lpdf with_cs_offsets
###   us_chol_cor expand_b aterms_for_newdata has_trunc as_priorlist
###   check_prior_slots resolve_prior_input neg_log_prior_fn resolve_bounds
###   spec_target spec_spelling frmtmb_register_prior_defaults ncp_eligible
###   ncp_scale_b ncp_unscale_b covstruct_has_chol block_sd_idx
###   block_cor_prior block_n_cor is_student_block sim_can sim_note
###   sim_context sim_draw sim_is_structured par_name_bare outer_par_names
###   estimated_coef_names log_sd_theta_index sdr_of require_fitted
###   hyp_parse_all hyp_vals_only hyp_env_vals hyp_eval hyp_tail_p
###   hyp_class_prefix hyp_labels hyp_samples_frame hyp_brms_result
###   brms_coef_names brms_re_rnames brms_block_has_r brms_coef_table
###   brms_stan_name brms_group_name brms_levels brms_re_parts hyp_eval_in
###   hyp_expr_vars brms_par_labels varcorr_matrices varcorr_layout
###   varcorr_values ce_grids_build ce_boot_one ce_frame ce_finalize
###   ce_cats_display ce_display_kind ce_pred_dpar ce_group_vars
###   ce_new_level_spec ce_boot_grids ce_draw_new_levels ce_structure_check
###   ce_re_formula ce_dots find_linpred arg_unset re_form_arg
###   frm_install_generics frm_check_dots

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# the objective seam: the bare likelihood closure of the frame
nll <- build_objective(fit$frame)
nll(fit$estimates)

# a block reader: this one has a diagonal factor and one sd, so it
# can be sampled non-centered
bk <- fit$frame$re_blocks[[1]]
c(eligible = ncp_eligible(bk), n_cor = block_n_cor(bk),
  cor_prior = block_cor_prior(bk))

# parameter labeling, in the two spellings
head(outer_par_names(fit))
head(par_name_bare(outer_par_names(fit)))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
