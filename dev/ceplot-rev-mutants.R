# Reviewer mutants (lane ceplot), injected into the loaded namespaces by
# dev/ceplot-rev-runtest.R before a test file runs. Each mutant edits the
# deparsed installed function and asserts its pattern matched exactly
# once, so a mutant that silently changed nothing cannot pass as caught.
mutants <- list(
  relabel_off = list(
    list("frmtmb", "ce_plan_eval",
         'fp$frame[["re_blocks"]][[id]][["levels"]] <- lv', "NULL")),
  relabel_wrongslot = list(
    list("frmtmb", "ce_plan_part",
         'idx = b$bk[["b_idx"]][seq_len(b$d)]',
         'idx = b$bk[["b_idx"]][b$d + seq_len(b$d)]')),
  relabel_badlabel = list(
    list("frmtmb", "ce_plan_part",
         'new[[i]]$label <- ce_design_label(blocks[[new[[i]]$block]], ',
         'new[[i]]$label <- "zz:zz"; list(')),
  valid_off = list(
    list("frmtmb", "ce_valid_effects", "if (all(ok))", "if (TRUE)")),
  valid_nogroup = list(
    list("frmtmb", "ce_valid_effects", ", cs_vars, ce_group_vars(x)))",
         ", cs_vars))")),
  valid_nocs = list(
    list("frmtmb", "ce_valid_effects", ", cs_vars, ce_group_vars(x)))",
         ", ce_group_vars(x)))")),
  valid_nodup = list(
    list("frmtmb", "ce_valid_effects", "&& !anyDuplicated(v)", "")),
  ol_off = list(
    list("frmtmb", "new_level_pick_apply", "j[hit] <- as.integer(p[hit])",
         "NULL")),
  ol_first = list(
    list("frmtmb", "old_level_pick_add",
         'sample.int(bl[[key]]$bk[["n_levels"]], \n            1L)', "1L")),
  ol_nofit = list(
    list("frmtmb", "fitted_old_levels",
         'object[["new_level_pick"]] <- pick', "NULL")),
  mmby_route_off = list(
    list("frmtmb", "ce_mm_in_block", "if (is.null(by)) ", "if (TRUE) ")),
  ordwarn_off = list(
    list("frmtmb", "ce_display_kind",
         'frm_warning("Predictions are treated as continuous variables in ",',
         'if (FALSE) frm_warning("Predictions are treated as continuous variables in ",')),
  mi_check_off = list(
    list("frmtmb", "ce_mi_idx_check",
         'for (lp in x$frame[["linpreds"]] %||% list()) {',
         'return(invisible(NULL))\n    for (lp in x$frame[["linpreds"]] %||% list()) {')),
  plot_false_draws = list(
    list("frmtmb", "plot.frmtmb_conditional_effects", "if (plot) {",
         "if (TRUE) {")),
  ps_order_off = list(
    list("frmtmb.sample", "ps_select", "if (length(bi) > 1L)",
         "if (FALSE)")),
  nsw_postwarmup = list(
    list("frmtmb.sample", "nsamples_saved",
         "as.integer(sim$n_save[1L] * sim$chains)",
         "as.integer((sim$n_save[1L] - sim$warmup2[1L]) * sim$chains)"))
)
apply_mutant <- function(name) {
  edits <- mutants[[name]]
  if (is.null(edits)) stop("unknown mutant ", name)
  for (ed in edits) {
    pkg <- ed[[1]]
    fn <- ed[[2]]
    ns <- asNamespace(pkg)
    f <- get(fn, envir = ns)
    txt <- paste(deparse(f), collapse = "\n")
    hits <- gregexpr(ed[[3]], txt, fixed = TRUE)[[1]]
    n <- sum(hits > 0)
    if (n != 1L) stop("mutant ", name, ": pattern matched ", n, " times in ",
                      fn)
    new <- eval(parse(text = sub(ed[[3]], ed[[4]], txt, fixed = TRUE)))
    environment(new) <- environment(f)
    targets <- list(ns)
    for (other in c("frmtmb.sample")) {
      if (isNamespaceLoaded(other)) {
        imp <- parent.env(asNamespace(other))
        if (exists(fn, envir = imp, inherits = FALSE)) {
          targets[[length(targets) + 1L]] <- imp
        }
      }
    }
    for (env in targets) {
      locked <- bindingIsLocked(fn, env)
      if (locked) unlockBinding(fn, env)
      assign(fn, new, envir = env)
      if (locked) lockBinding(fn, env)
    }
    # an S3 method registered for dispatch is looked up in the S3 table
    if (grepl("^(plot|print)[.]", fn)) {
      cls <- sub("^(plot|print)[.]", "", fn)
      gen <- sub("[.].*$", "", fn)
      registerS3method(gen, cls, new, envir = ns)
    }
  }
  invisible(TRUE)
}
