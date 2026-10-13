# Reviewer of the 0.70.0 consolidation: sourced by a study script after
# library(frmtmb) when RELREV_INSTR=1. Wraps obj_settle() so that each
# call records who called it, whether the kept optimum is the start, and
# whether restoring the saved state alone (lane moparity's restore, the
# code the merge dropped) would already leave the objective at the kept
# point. The wrapper only reads; the original then runs unchanged.
relrev_log <- new.env()
relrev_log$rows <- list()
local({
  ns <- asNamespace("frmtmb")
  orig <- get("obj_settle", ns)
  w <- function(obj, opt, best, saved) {
    caller <- tryCatch(deparse(sys.call(-1L)[[1L]])[1L],
                       error = function(e) "?")
    same <- identical(best, opt)
    at_restore <- NA
    if (same) {
      get("obj_state_restore", ns)(obj, saved)
      at_restore <- get("obj_state_at", ns)(obj, best$par)
    }
    r <- orig(obj, opt, best, saved)
    at_end <- get("obj_state_at", ns)(obj, best$par)
    relrev_log$rows[[length(relrev_log$rows) + 1L]] <- data.frame(
      caller = caller, same = same, at_restore = at_restore,
      at_end = at_end, seed = relrev_log$seed %||% NA_integer_)
    r
  }
  environment(w) <- environment()
  utils::assignInNamespace("obj_settle", w, "frmtmb")
})
relrev_summary <- function(tag) {
  d <- do.call(rbind, relrev_log$rows)
  if (is.null(d)) {
    cat("INSTR", tag, ": no obj_settle() calls\n")
    return(invisible(NULL))
  }
  tb <- as.data.frame(table(caller = d$caller, same = d$same,
                            at_restore = d$at_restore, at_end = d$at_end,
                            useNA = "ifany"))
  tb <- tb[tb$Freq > 0, ]
  cat("INSTR", tag, ": obj_settle() calls", nrow(d), "\n")
  print(tb, row.names = FALSE)
  saveRDS(d, file.path("dev/relrev070-log", paste0("instr-", tag, ".rds")))
}
