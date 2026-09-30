# Reviewer mutants for the punch-round-1 re-check (lane ceplot): my 16
# from dev/ceplot-rev-mutants.R, four re-pointed at the functions the
# punch round rewrote, plus the worker's 11 from dev/ceplot-p1-mutants.R
# under a "w_" prefix. apply_mutant() is the one of the first file.
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-mutants.R")
mine <- mutants
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-p1-mutants.R")
theirs <- mutants
mine$valid_nogroup <- list(
  list("frmtmb", "ce_model_vars", "unique(c(v, by, ce_group_vars(x)))",
       "unique(c(v, by))"))
mine$valid_nocs <- list(
  list("frmtmb", "ce_model_vars",
       'unlist(lapply(lp[["cs"]] %||% list(), function(ct) {',
       "unlist(lapply(list(), function(ct) {"))
mine$ol_off <- list(
  list("frmtmb", "new_level_pick_apply", "j[hit] <- jj[hit]", "NULL"))
mine$ol_first <- list(
  list("frmtmb", "old_level_pick_add", "lv[sample.int(length(lv), 1L)]",
       "lv[1L]"))
mutants <- c(mine, stats::setNames(theirs, paste0("w_", names(theirs))))
