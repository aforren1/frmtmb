# Shared harness for the hand-translated brms vignettes.
#
# ATTACH frmtmb.sample BEFORE RUNNING THESE SCRIPTS. Several of them
# call frm_sample(), check_laplace(), log_lik() or loo(), which moved to
# extensions/frmtmb.sample when the sampling surface was split out
# (dev/draws-extraction.md). The scripts are deliberately NOT rewritten
# to qualify those names: they are hand translations of brms vignettes
# and the point is that the translated line reads like the brms line.
# Install the extension and `library(frmtmb.sample)` first, and they run
# as written.
#
# Why a harness and not plain `try()`: the audit needs a machine-countable
# record of what each call did, so the scoreboard in
# `dev/brms-vignette-audit.md` is measured and not remembered. Every
# translated call goes through `bv()`, which runs it, catches the error,
# and appends one row to a CSV the summarizer reads.
#
# Edge vocabulary (the `edge` argument), applied to the TRANSLATION, not
# to the outcome:
#   NA         the brms line runs unchanged apart from brm -> frm and the
#              removal of MCMC-only arguments
#   "SPELLING" works, under a different name or argument
#   "BEHAVIOR" runs, but the output differs from what brms documents
#   "MISSING"  no frmtmb path
#   "REFUSAL"  frmtmb refuses on purpose

bv_state <- new.env(parent = emptyenv())

bv_init <- function(vignette) {
  bv_state$vignette <- vignette
  bv_state$rows <- list()
  bv_state$t0 <- Sys.time()
  out <- Sys.getenv("BV_OUT", unset = tempdir())
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  bv_state$csv <- file.path(out, paste0(vignette, ".csv"))
  cat("### vignette:", vignette, "\n")
  invisible(NULL)
}

#' Run one translated call.
#'
#' @param kind "model" for a frm()/frm_multiple() call, "post" for
#'   anything downstream, "data" for setup that is not part of the claim.
#' @param label the vignette's own name for the call.
#' @param expr the translated expression.
#' @param edge one of the vocabulary above, or NA for a clean port.
#' @param why one line saying what the divergence is.
#' @return the value, or NULL when the call failed.
bv <- function(kind, label, expr, edge = NA_character_, why = "") {
  t0 <- Sys.time()
  cat("\n--- [", kind, "] ", label,
      if (!is.na(edge)) paste0("  <", edge, "> ", why) else "", "\n", sep = "")
  val <- tryCatch(expr,
                  error = function(e) structure(conditionMessage(e),
                                                class = "bv_error"))
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
  ok <- !inherits(val, "bv_error")
  if (!ok) cat("ERROR: ", as.character(val), "\n", sep = "")
  # the log is the evidence for the BEHAVIOR edges, so a post-processing
  # call's value is printed. The null device keeps a print method that
  # draws (conditional_effects, ggplot) from needing a real one.
  if (ok && kind != "model" && !is.null(val) && !isTRUE(bv_state$quiet)) {
    grDevices::pdf(NULL)
    try(print(val), silent = TRUE)
    grDevices::dev.off()
  }
  cat("[", if (ok) "ok" else "err", " ", secs, "s]\n", sep = "")
  bv_state$rows[[length(bv_state$rows) + 1L]] <- data.frame(
    vignette = bv_state$vignette, kind = kind, label = label,
    edge = edge, why = why, ok = ok, secs = secs,
    msg = if (ok) "" else as.character(val),
    stringsAsFactors = FALSE)
  if (ok) invisible(val) else invisible(NULL)
}

bv_done <- function() {
  d <- do.call(rbind, bv_state$rows)
  utils::write.csv(d, bv_state$csv, row.names = FALSE)
  tot <- round(as.numeric(difftime(Sys.time(), bv_state$t0, units = "secs")), 1)
  cat("\n### done:", bv_state$vignette, nrow(d), "calls,", tot, "s\n")
  cat("### models:", sum(d$kind == "model"), " ok:",
      sum(d$kind == "model" & d$ok), "\n")
  cat("### post:", sum(d$kind == "post"), " ok:",
      sum(d$kind == "post" & d$ok), "\n")
  print(table(d$kind, ifelse(is.na(d$edge), "CLEAN", d$edge)))
  invisible(d)
}

# frmtmb is loaded from the worktree, brms is NOT attached: a
# post-processing name must resolve to frmtmb or fail, which is part of
# what is being measured. brms data and reference values are reached
# with the `brms::` prefix. `brms-coexistence.R` is the one script that
# attaches brms on purpose.
bv_load <- function() {
  suppressMessages(pkgload::load_all(
    Sys.getenv("BV_PKG", unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig"),
    quiet = TRUE, export_all = FALSE))
  options(warn = 1)
}
