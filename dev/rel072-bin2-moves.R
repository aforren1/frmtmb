# The verdict moves of the 0.72.0 consolidation in bin 2 of the ported
# brms suite. dev/brmsport-bin2-verdicts-2p.R (lane bin2port's) sources
# this just before it writes its two manual tables, as the 0.71.0
# consolidation's dev/rel071-bin2-moves.R did for its rows; the
# resolve (dev/rel072-resolve.R) adds that call. Each move is a row the
# merged build records differently from the lane's, with the reason.
# First run on the merged build: dev/rel072-log/bin2-1.txt, the ledger
# stopping on exactly these 7 rows.
#
# Returns function(X, kind): kind "2p" edits the row verdicts, kind
# "values" the value-check verdicts.
function(X, kind) {
  drop <- function(X, ids) {
    stopifnot(all(ids %in% X$id))
    cat("bin-2", kind, "verdicts removed (now hold):", ids, "\n")
    X[!X$id %in% ids, , drop = FALSE]
  }
  set <- function(X, id, verdict, class, reason, msg_has = NULL) {
    i <- which(X$id == id)
    stopifnot(length(i) == 1L)
    X$verdict[i] <- verdict
    X$class[i] <- class
    X$reason[i] <- reason
    if (!is.null(msg_has)) X$msg_has[i] <- msg_has
    cat("bin-2", kind, "verdict moved:", id, "to", verdict, "/", class, "\n")
    X
  }
  if (kind == "2p") {
    # lane catmulti's bare multinomial() reads its categories from the
    # response, so the adapter's family (the bare constructor, resolved
    # on brms's ncat categories) now exists and both rows hold;
    # extraDistr is installed, so brms's reference dbeta_binomial() runs
    # and log_lik:267 holds
    X <- drop(X, c("log_lik:504", "posterior_predict:349", "log_lik:267"))
    # (two rows set to "absent" here, log_lik:510 and
    # posterior_predict:355, hold since 0.73.0: dev/rel073-bin2-moves.R)
  } else {
    # extraDistr is installed: the three references run, and hold
    X <- drop(X, c("distributions:123", "distributions:177",
                   "log_lik:436"))
    X <- set(X, "log_lik:462", "divergence", "domain", paste(
      "cumulative with all three thresholds 0 (brms's prep): after lane",
      "catmulti's fix of ord_cumulative_logpmf() frmtmb's end",
      "categories are log F(-disc mu) and log(1 - F(-disc mu)) to",
      "4.4e-16, and the two zero-width middle categories are -Inf, the",
      "exact value, where brms's log_lik_cumulative() gives NaN",
      "(dev/rel072-cumtie-probe.R: brms NaN on 200 of 400 values, the",
      "zero-width cells; frmtmb NaN on none). As log_lik:478"),
      msg_has = "NA in 0 values from the adapter, 200 from brms")
  }
  X
}
