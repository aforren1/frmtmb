# The verdict moves of the 0.81.0 consolidation in bin 2T of the ported
# brms suite. dev/brmsport-bin2-verdicts-2t.R sources this after the
# 0.71.0 and 0.73.0 moves.
#
# Why a file: lane hollow edited dev/brmsport-bin2-verdicts-manual.tsv
# by hand (its Row 2), but dev/brmsport-bin2-all.sh rewrites that file
# from -2t.R, so the consolidation's regeneration put the 0.80.0 rows
# back and brmsfit-methods:469 stayed divergence (hollow). These are
# hollow's three rows exactly as its tree has them
# (dev/hollow-findings.md, "Row 2" and the fix round).
#
# Returns function(X, kind): kind "2t" edits the 2T table.
function(X, kind) {
  stopifnot(identical(kind, "2t"))
  S <- "frmtmb.sample"
  m <- function(k) paste0("brmsfit-methods:", k)
  # draws carry $data, so loo(fit6, newdata = fit6$data) is the elpd on
  # the fitted rows, which is brms's: the row holds and has no verdict
  sel <- X$id == m(469) & X$pkg == S
  stopifnot(sum(sel) == 1L)
  X <- X[!sel, , drop = FALSE]
  cat("bin-2 2t verdicts removed (now hold):", m(469), S, "\n")
  for (k in c(508, 544)) {
    i <- which(X$id == m(k) & X$pkg == S)
    stopifnot(length(i) == 1L)
    X$reason[i] <- paste(
      "the loo_*() functions refuse newdata other than the fitted rows, as",
      "loo() does: a leave-one-out weight is the observation's own column",
      "of log_lik(), so it exists only for a row the model was fitted on.",
      "brms computes the weights on the new rows' own response (this row",
      "with allow_new_levels = TRUE). Lane loodraws's decision 1")
    cat("bin-2 2t reason edited:", m(k), S, "\n")
  }
  X
}
