# The verdict moves of the 0.71.0 consolidation to bin 2 of the ported
# brms suite (lane remeasure's port). Lane remeasure's verdict scripts,
# dev/brmsport-bin2-verdicts-2t.R and -2p.R, write the manual tables
# from R; dev/rel071-resolve.R makes each pass its table through the
# function below before writing it, so the moves survive every rerun of
# dev/brmsport-bin2-all.sh. The bin-2 ledger stops on a verdict for an
# assertion that holds and on a "hollow" row that no longer holds; on
# the merged build it named the 14 rows of the first two groups
# (dev/rel071-log/bin2-six.txt, bin2.txt). The last two groups are
# reasons that the round made false.
#
# 2T, lane postfix (item h of dev/postfix-findings.md):
# - eight rows of frmtmb.sample's methods file hold (bayes_R2(newdata =),
#   logLik() of draws, loo() and waic() of several models);
# - brmsfit-methods:469 held hollowly on a NULL newdata that loo()
#   swallowed, and now fails at its setup on the refusal of `newdata`
#   by name: a divergence, policy.
# 2P, lane famuni:
# - log_lik:271, :324 and posterior_predict:159, :184, :190 hold through
#   the existing adapter now that discrete_weibull(), frechet() and
#   gen_extreme_value() exist (log_lik:324 compares the frechet log
#   density with brms's dfrechet() value by value; the other four check
#   lengths, as brms wrote them);
# - log_lik:349 fails on rows outside the GEV support, where frmtmb's
#   log density is finite and steep and brms's dgen_extreme_value() is
#   NaN; inside the support they agree (dev/rel071-gev-probe.R,
#   dev/rel071-log/gev-probe.txt): a divergence, domain.
# Reasons the round made false (each row still fails; the class says
# why):
# - the dirichlet, dirichlet2 and logistic_normal rows of the categorical
#   blocks: lane famsimplex added these families; the thin adapter does
#   not map brms's per-category mu or the logistic-normal lncor, as for
#   categorical (records: "prep$dpars has no 'mu'", "the prep carries
#   lncor"): not ported (cost). The dirichlet_multinomial rows stay
#   absent;
# - the G1 frechet and gen_extreme_value rows (distributions:92 to 141):
#   lane famuni added the families, but brms's d/p/r functions of them
#   are not shimmed (`brms_g1_shims` covers the families it was built
#   for): not ported (cost).
function(V, part) {
  set_row <- function(V, id, pkg, verdict, class, reason) {
    i <- which(V$id == id & V$pkg == pkg)
    stopifnot(length(i) == 1L)
    V$verdict[i] <- verdict
    V$class[i] <- class
    V$reason[i] <- reason
    V
  }
  if (part == "2t") {
    S <- "frmtmb.sample"
    hold <- paste0("brmsfit-methods:", c(121, 434, 446, 447, 448, 449,
                                         1015, 1016))
    stopifnot(all(hold %in% V$id[V$pkg == S]))
    V <- V[!(V$id %in% hold & V$pkg == S), , drop = FALSE]
    V <- set_row(V, "brmsfit-methods:469", S, "divergence", "policy", paste(
      "loo() of draws refuses newdata by name, as log_lik() of draws",
      "does: draws are evaluated on the fitted objective's frame, and a",
      "frame rebuilt from newdata would rebuild a smooth's basis on it",
      "(lane postfix, item h). brms computes the elpd on newdata. The",
      "row held hollowly at 0.70.0, on a NULL newdata that loo()",
      "swallowed"))
  }
  if (part == "2p") {
    F <- "frmtmb"
    hold <- c("log_lik:271", "log_lik:324", "posterior_predict:159",
              "posterior_predict:184", "posterior_predict:190")
    stopifnot(all(hold %in% V$id))
    V <- V[!V$id %in% hold, , drop = FALSE]
    V <- set_row(V, "log_lik:349", F, "divergence", "domain", paste(
      "frmtmb's gen_extreme_value() log density is finite and steep",
      "outside the support (-1.0e150 at y = 0, mu 1, sigma 0.1, xi 0.3),",
      "so the optimizer is turned back and a mixture component adds",
      "nothing there; brms's dgen_extreme_value() is NaN. Inside the",
      "support the two agree (-1.55396; dev/rel071-gev-probe.R). Lane",
      "famuni's documented choice"))
    for (id in c("log_lik:516", "log_lik:523", "log_lik:534",
                 "posterior_predict:360", "posterior_predict:361",
                 "posterior_predict:368", "posterior_predict:369",
                 "posterior_predict:380", "posterior_predict:381")) {
      V <- set_row(V, id, F, "cannot transfer", "not-ported", paste(
        "dirichlet, dirichlet2 and logistic_normal are frmtmb families",
        "since 0.71.0 (lane famsimplex), but brms's prep holds one",
        "predictor per non-reference category (and logistic_normal's",
        "lncor), which the thin adapter does not map ('prep$dpars has no",
        "mu', 'the prep carries lncor'), as for categorical: transferable",
        "at a cost"))
    }
    for (id in c("log_lik:510", "posterior_predict:355")) {
      V <- set_row(V, id, F, "cannot transfer", "absent", paste(
        "dirichlet_multinomial is not a frmtmb family ('not a supported",
        "family'): a gap, not a decision"))
    }
    for (k in c(92, 95, 97, 100, 102, 134, 137, 139, 141)) {
      V <- set_row(V, paste0("distributions:", k), F, "cannot transfer",
        "not-ported", paste(
        "frechet and gen_extreme_value are frmtmb families since 0.71.0",
        "(lane famuni), but brms's d/p/r functions of them are not",
        "shimmed (brms_g1_shims covers the families it was built for), so",
        "the row stops at 'could not find function': transferable at a",
        "cost"))
    }
  }
  V
}
