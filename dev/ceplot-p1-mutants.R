# Lane ceplot punch 1: mutants of this round's guards, injected into the
# loaded namespaces before a test file runs (the construction of
# dev/ceplot-rev-mutants.R, whose apply_mutant() is reused). Each edits
# the deparsed installed function and asserts its pattern matched
# exactly once, so a mutant that changed nothing cannot pass as caught.
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-mutants.R")
mutants <- list(
  relabel_wrongslot = list(
    list("frmtmb", "ce_plan_part",
         'idx = b$bk[["b_idx"]][seq_len(b$d)]',
         'idx = b$bk[["b_idx"]][b$d + seq_len(b$d)]')),
  rename_na_off = list(
    list("frmtmb", "ce_plan_part", "if (anyNA(vapply(b$gv",
         "if (FALSE && anyNA(vapply(b$gv")),
  valid_nocs = list(
    list("frmtmb", "ce_model_vars",
         'all.vars(ct[["expr"]] %||% str2lang(sub("^cs", "", ',
         "c(c(c(")),
  valid_nlpars = list(
    list("frmtmb", "ce_model_vars", "c(r$nlpars, names(r$dpars))",
         "character(0)")),
  ol_perblock = list(
    list("frmtmb", "new_level_key", 'bk[["group_name"]] %||% ',
         "NULL %||% ")),
  ol_offrows = list(
    list("frmtmb", "old_level_pick_add", "which(rp$is_new %||% is.na(rp$j))",
         "which(is.na(rp$j))")),
  group_values_off = list(
    list("frmtmb", "group_values", 'identical(expr[[1L]], as.name(":"))',
         "FALSE")),
  rawvars_off = list(
    list("frmtmb", "ce_base_frame", "if (is.null(raw) ||", "if (TRUE ||")),
  predict_after_warn = list(
    list("frmtmb", "ce_display_kind", 'if (identical(method, "predict"))',
         "if (FALSE)")),
  plot_ignored_off = list(
    list("frmtmb", "ce_plot_ignored", "if (length(nms))", "if (FALSE)")),
  ps_nl_off = list(
    list("frmtmb.sample", "ps_select", "!is_nl", "TRUE"))
)
