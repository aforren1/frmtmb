# Review of lane adhess, item 6: frm_coherence() near the edges, where
# the refusal's condition is absent (an interval inside (0, 1)) and
# where it is present; and the mirror image near a coherence of 0.
#   REV_ARM=lane|base Rscript dev/adhess-rev-coh.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
suppressMessages(library(frmtmb.coupling))
cat("coupling from", dirname(find.package("frmtmb.coupling")), "\n")
tf <- paste0("C:/Users/adf44/source/r/frmtmb-wt-adhess/extensions/",
             "frmtmb.coupling/tests/testthat/test-coherence.R")
ex <- parse(tf)
for (e in ex[1:3]) eval(e)
for (eta in c(-16, -8, 3, 6, 9, 12, 16)) for (nseg in c(4L, 32L)) {
  set.seed(301)
  d <- cp_units(rep(eta, 15L), rep(0.4, 15L), n = nseg)
  r <- cap(frm(cp_bf(), family = cross_wishart(), data = d))
  fit <- r$v
  if (inherits(fit, "rev_err")) {
    cat(sprintf("eta %3g n %2d | fit ERROR %s\n", eta, nseg, fit))
    next
  }
  co <- tryCatch(frm_coherence(fit, newdata = d[1, ]),
                 error = function(e) conditionMessage(e))
  if (is.character(co)) {
    cat(sprintf("eta %3g n %2d | code %d | REFUSED: %s\n", eta, nseg,
                fit$opt$convergence, substr(co, 1, 90)))
  } else {
    cat(sprintf(paste0("eta %3g n %2d | code %d | est %.10g se %.3g | ",
                       "[%.6g, 1 - %.3g] | warn %d\n"),
                eta, nseg, fit$opt$convergence, co$.eta, co$.se, co$.lower,
                1 - co$.upper, length(r$w)))
  }
}
