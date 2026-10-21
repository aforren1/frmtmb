# Reviewer, fix round 2: the gaussian [9, 9.5] KS check of cdfs-rev-r2.R
# with the truncated CDF written correctly (the first version negated it)
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(3)
famg <- ns$as_frmtmb_family(gaussian())
for (w in list(c(9, 9.5), c(30, 31), c(-0.5, 0.5))) {
  sv <- ns$sim_response(famg, list(mu = rep(0, 4000), sigma = rep(1, 4000)),
                        list(trunc_lb = rep(w[1], 4000), trunc_ub = rep(w[2], 4000)), 4000)
  S <- function(q) pnorm(q, lower.tail = FALSE, log.p = TRUE)
  Ft <- function(q) (1 - exp(S(q) - S(w[1]))) / (1 - exp(S(w[2]) - S(w[1])))
  cat(sprintf("gaussian in [%g, %g]: outside %d, KS p %.3g\n", w[1], w[2], sum(sv < w[1] | sv > w[2]),
              suppressWarnings(ks.test(sv, Ft)$p.value)))
}
