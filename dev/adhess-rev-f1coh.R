# Re-check of fix round 1, finding 4: the symmetric logit-scale bound,
# by mocking cp_link_se() at both ends, inside and past the bound.
#   REV_ARM=lane Rscript dev/adhess-rev-f1coh.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
suppressMessages(library(frmtmb.coupling))
cns <- asNamespace("frmtmb.coupling")
edge <- cns$cp_logit_edge
cat(sprintf("edge %.6f, plogis(edge) = 1 - %.3g, plogis(-edge) = %.3g\n",
            edge, 1 - plogis(edge), plogis(-edge)))
z <- stats::qnorm(0.975)
for (side in c(-1, 1)) for (off in c(-0.1, 0.1)) {
  eta <- side * 30
  se <- (edge - 30) / z + off
  r <- tryCatch(testthat::with_mocked_bindings(
    frm_coherence(NULL),
    cp_link_se = function(...) list(fit = eta, se.fit = se, z = z),
    .package = "frmtmb.coupling"), error = function(e) conditionMessage(e))
  cat(sprintf("eta %+g se %.4f (%s the bound by 0.1 SE) | ", eta, se,
              if (off > 0) "past" else "inside"))
  if (is.character(r)) cat("REFUSED\n") else
    cat(sprintf("[%.3g, 1 - %.3g] inside (0,1): %s\n", r$.lower,
                1 - r$.upper, r$.lower > 0 && r$.upper < 1))
}
