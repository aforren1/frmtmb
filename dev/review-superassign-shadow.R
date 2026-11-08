# Review of lane superassign: inside each changed function, find every
# binding (a plain assignment, a for variable, a formal of a nested
# function) of the name the lane gave its new environment, other than
# the defining `name <- new.env(...)`. A nested binding of the same name
# would make `name$field` read the wrong object in that scope.
lane_root <- "C:/Users/adf44/source/r/frmtmb-wt-superassign"
envs <- list(
  "R/autoscale.R" = c(autoscale_choose = "got"),
  "R/brms-shapes.R" = c(fit_fd_se = "geta"),
  "R/ce-levels.R" = c(ce_level_plan = "tally"),
  "R/compat.R" = c(compat_rule_builder = "acc", compat_hand_rules_tbl = "acc"),
  "R/confint.R" = c(par_alias_index = "got", varcorr_trans_rows = "acc",
                    hyp_env_vals = "acc"),
  "R/curve-cov.R" = c(curve_catch_span = "got"),
  "R/fit.R" = c(nlminb_trial_fn = "st", newton_polish = "rp",
                mo_search = "tried"),
  "R/frame.R" = c(resolve_group_calls = "acc", frame_raw_vars = "seen",
                  assemble_frame = "rhs_acc", assemble_frame = "codes_acc"),
  "R/generic-owners.R" = c(frm_bind_generic = "memo"),
  "R/gp-by.R" = c(gp_brms_natural = "acc"),
  "R/hazard-containers.R" = c(frm_hazard_reads = "acc"),
  "R/importance.R" = c(build_importance_objective = "shift"),
  "R/me.R" = c(me_newdata_extra = "obs", me_node_value = "memo"),
  "R/methods-fit.R" = c(VarCorr.frmtmb_fit = "cur"),
  "R/par-template.R" = c(nl_restart_round = "h1", nl_fit_off_saddle = "hq"),
  "R/parse.R" = c(mo_bar_rewrite = "acc", expand_double_verts = "queue",
                  parse_linpred = "track", rhs_split_ac = "found"),
  "R/priors.R" = c(as_priorlist = "acc", prior_table = "acc",
                   cov_class_coefs = "acc", resolve_priorlist = "st"),
  "R/ps.R" = c(ps_extract = "found"),
  "R/re-formula.R" = c(re_formula_bars = "found"),
  "R/simulate-new.R" = c(nat_slots = "acc"),
  "extensions/frmtmb.ode/R/ode.R" = c(frm_ode = "got"),
  "extensions/frmtmb.sample/R/loo.R" = c(log_lik.frmtmb_draws = "seen"),
  "extensions/frmtmb.sample/R/methods-draws.R" = c(draws_laplace_watch = "seen"),
  "extensions/frmtmb.sample/R/prior-draws.R" = c(prior_plan = "acc"),
  "extensions/frmtmb.sample/R/sample.R" = c(default_priors_for = "acc"))

bindings <- function(e, name, depth = 0L, acc = list()) {
  if (is.call(e)) {
    hd <- e[[1L]]
    if (identical(hd, as.name("function"))) {
      if (name %in% names(e[[2L]])) {
        acc[[length(acc) + 1L]] <- paste0("formal of nested function (depth ",
                                          depth + 1L, ")")
      }
      return(bindings(e[[3L]], name, depth + 1L, acc))
    }
    if (as.character(hd)[1L] %in% c("<-", "=", "<<-") &&
          identical(e[[2L]], as.name(name))) {
      rhs <- deparse1(e[[3L]])
      if (!grepl("^new.env\\(parent = emptyenv\\(\\)\\)$", rhs) || depth > 0L) {
        acc[[length(acc) + 1L]] <- paste0("assign at depth ", depth, ": ",
                                          substr(rhs, 1, 60))
      }
    }
    if (identical(hd, as.name("for")) && identical(e[[2L]], as.name(name))) {
      acc[[length(acc) + 1L]] <- paste0("for variable at depth ", depth)
    }
    for (a in as.list(e)[-1L]) {
      if (!missing(a)) acc <- bindings(a, name, depth, acc)
    }
  }
  acc
}
cat("lib: none (source scan)\n")
n_hit <- 0L
for (fl in names(envs)) {
  ex <- parse(file.path(lane_root, fl), keep.source = FALSE)
  tops <- list()
  for (e in ex) {
    if (is.call(e) && as.character(e[[1L]])[1L] %in% c("<-", "=") &&
          is.name(e[[2L]])) {
      tops[[as.character(e[[2L]])]] <- e[[3L]]
    }
  }
  sel <- envs[[fl]]
  for (i in seq_along(sel)) {
    fn <- names(sel)[i]
    nm <- sel[[i]]
    b <- bindings(tops[[fn]], nm)
    # the defining `name <- new.env()` sits at depth 1 (inside the
    # function); a second one at depth 1 is the case to read
    defs <- b[grepl("^assign at depth 1: new.env", unlist(b))]
    other <- setdiff(unlist(b), unlist(defs))
    cat(sprintf("%-28s %-10s definitions %d, other bindings %d\n", fn, nm,
                length(defs), length(other)))
    for (o in other) cat("   ", o, "\n")
    n_hit <- n_hit + length(other)
  }
}
cat("RESULT other bindings of the new names:", n_hit, "\n")
