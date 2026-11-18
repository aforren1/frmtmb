# rel084-kt copy of dev/slopetrig-grby.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed); the before arm tsv is read from the slopetrig lane tree by absolute path. No other change.
# Lane slopetrig, fix round: the review's gr(by =) year design
# (dev/reviews/slopetrig-review-designs.R, grby_year), seeds 1 to 10,
# after the ratio is read over the rows the block covers. Arms in one
# process: ref (centered, standardized column), F (autoscale = FALSE),
# T (autoscale = TRUE), D (default). Prints per seed whether D engaged,
# each arm's shortfall against the best arm and D's warning count, then
# the counts. The before arm is the review's own run of the previous
# lane build (dev/reviews/slopetrig-review-oos-lane.tsv).
# Usage: Rscript slopetrig-grby.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb lib:", dirname(find.package("frmtmb")), "\n")
source("dev/reviews/slopetrig-review-designs.R")
ds <- designs$grby_year
run1 <- function(expr) {
  nw <- 0L
  f <- withCallingHandlers(expr,
    warning = function(w) {
      nw <<- nw + 1L
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  list(ll = as.numeric(logLik(f)), nw = nw, as = isTRUE(f$autoscaled))
}
rows <- list()
for (seed in 1:10) {
  d <- ds$data(seed)
  ref <- run1(frm(ds$form, data = ds$ref(d), family = ds$family))
  rF <- run1(frm(ds$form, data = d, family = ds$family,
                 control = frmtmb_control(autoscale = FALSE)))
  rT <- run1(frm(ds$form, data = d, family = ds$family,
                 control = frmtmb_control(autoscale = TRUE)))
  rD <- run1(frm(ds$form, data = d, family = ds$family))
  best <- max(ref$ll, rF$ll, rT$ll, rD$ll)
  rows[[seed]] <- data.frame(seed = seed, engaged = rD$as,
                             sh_ref = signif(best - ref$ll, 3),
                             sh_F = signif(best - rF$ll, 3),
                             sh_T = signif(best - rT$ll, 3),
                             sh_D = signif(best - rD$ll, 3), w_D = rD$nw)
}
r <- do.call(rbind, rows)
print(r, row.names = FALSE)
cat(sprintf("after: engaged %d of 10 ; D short %d (silent %d, max %.3g) ; T short %d ; ref short %d\n",
            sum(r$engaged), sum(r$sh_D > 1e-3),
            sum(r$sh_D > 1e-3 & r$w_D == 0), max(r$sh_D),
            sum(r$sh_T > 1e-3), sum(r$sh_ref > 1e-3)))
b <- utils::read.delim("C:/Users/adf44/source/r/frmtmb-wt-slopetrig/dev/reviews/slopetrig-review-oos-lane.tsv")
b <- b[b$design == "grby_year", ]
bb <- pmax(b$ll_ref, b$ll_F, b$ll_T, b$ll_D)
cat(sprintf("before (review's run of the previous build): engaged %d of %d ; D short %d (silent %d, max %.3g)\n",
            sum(b$as_D), nrow(b), sum(bb - b$ll_D > 1e-3),
            sum(bb - b$ll_D > 1e-3 & b$w_D == 0), max(bb - b$ll_D)))
