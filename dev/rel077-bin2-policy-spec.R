# The rows that test brms's own numeric helper API: a divergence by
# policy since the 0.77.0 consolidation, where they were "cannot
# transfer / decision" (classes no-density and no-quantile). Six run
# through the adapter (dev/brmsport-bin2-verdicts-2p.R sets their
# verdict); the twelve of the multivariate normal and t block are not
# ported, and dev/brmsport-bin2-ledger.R gives them this verdict in
# place of the survey's class.
bin2_numeric_helper_reason <- paste(
  "divergence by policy: the row tests brms's own numeric helper API",
  "(dmulti_normal() and dmulti_student_t() with check =, rmulti_normal(),",
  "rmulti_student_t(), qstudent_t(), qskew_normal(), qfrechet(),",
  "qgen_extreme_value()).",
  "brms itself, mvtnorm and extraDistr provide them; frmtmb does not",
  "re-export them (SPEC.md: numeric d/p/q/r functions come from other",
  "packages)")
bin2_policy_unported <- paste0("distributions:",
  c(15, 17, 19, 20, 23, 25, 28, 31, 33, 36, 39, 42))
# A literal pass that a reader could take for more than it shows. The
# two rescor posterior_predict rows are frmtmb's (core) rows, as every
# bin-2P row is; frmtmb.sample's posterior_predict() of set_rescor(TRUE)
# draws is another path, with a defect filed at the 0.77.0 release.
bin2_pass_note <- stats::setNames(rep(paste(
  "pass through core's predict() joint draw",
  "(brms_adapt_mv_posterior_predict()), the route of every bin-2P row;",
  "it says nothing about",
  "frmtmb.sample's posterior_predict() of set_rescor(TRUE) draws, which",
  "draws the responses independently (dev/test-backlog.md, filed at the",
  "0.77.0 release)"), 2L),
  c("posterior_predict:69", "posterior_predict:72"))
