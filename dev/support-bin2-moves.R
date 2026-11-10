# The verdict moves of lane support (round of 2026-11-09) in bin 2 of the
# ported brms suite, in the form of dev/rel072-bin2-moves.R.
# dev/brmsport-bin2-verdicts-2p.R sources this just before it writes its
# two manual tables, so that a regeneration keeps the rows this lane
# edited in dev/brmsport-bin2-verdicts-manual-2p.tsv and
# dev/brmsport-bin2-valuechecks-manual.tsv. dev/support-findings.md has
# the measurements.
#
# Returns function(X, kind): kind "2p" edits the row verdicts, kind
# "values" the value-check verdicts.
function(X, kind) {
  if (kind == "2p") {
    # the GEV's log density off the tape is exact: -Inf outside the
    # support, where brms's is NaN
    i <- which(X$id == "log_lik:349")
    stopifnot(length(i) == 1L)
    X$verdict[i] <- "divergence"
    X$class[i] <- "domain"
    X$reason[i] <- paste(
      "frmtmb's gen_extreme_value() log density off the tape, which",
      "log_lik() reads, is -Inf outside the support, the exact value;",
      "brms's dgen_extreme_value() is NaN there. At the block's seed 3 of",
      "25 draws put y = 0 outside the support: frmtmb -Inf on all 3, brms",
      "NaN on all 3; on the 22 inside the two agree to 2.7e-15 relative",
      "(dev/support-gev349-probe.R). The fitted objective keeps its",
      "finite, steep barrier outside the support, which turns the",
      "optimizer back. As log_lik:462 and :478")
    cat("bin-2", kind, "verdict moved: log_lik:349 to divergence / domain\n")
  } else {
    # betareg 3.2-6 is installed: brms's log_lik_xbeta() runs, and the
    # check holds (50 values, max rel 6.4e-14)
    stopifnot("log_lik:580" %in% X$id)
    X <- X[X$id != "log_lik:580", , drop = FALSE]
    cat("bin-2", kind, "verdicts removed (now hold): log_lik:580\n")
  }
  X
}
