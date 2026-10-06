# Reviewer of lane ordmix, final check: section (3) of dev/ordmix-rev3-misc.R
# with the Hessian eigenvalues printed, by optimHess and by the tape.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib", "C:/Users/adf44/source/r/rellib-r5", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== (3) B4: identified count against the Hessian rank, no predictor\n")
rank_of <- function(f) {
  H <- optimHess(f$opt$par, f$obj$fn, f$obj$gr)
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  evt <- eigen(f$obj$he(f$opt$par), symmetric = TRUE, only.values = TRUE)$values
  cat("   eig/max (optimHess):", format(signif(sort(abs(ev), TRUE) / max(abs(ev)), 2)), "
")
  cat("   eig/max (tape he):  ", format(signif(sort(abs(evt), TRUE) / max(abs(evt)), 2)), "
")
  cat("   max abs gradient", signif(max(abs(f$obj$gr(f$opt$par))), 3), "
")
  sum(ev > 1e-6 * max(ev))
}
set.seed(20261007)
n <- 1500
g <- factor(sample(c("a", "b"), n, TRUE))
y4 <- sample(1:4, n, TRUE, prob = c(0.3, 0.2, 0.25, 0.25))
y3 <- sample(1:3, n, TRUE, prob = c(0.4, 0.35, 0.25))
d <- data.frame(g = g, y = ifelse(g == "a", y4, y3))
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
d$yh2 <- ifelse(runif(n) < ifelse(d$g == "a", 0.1, 0.35), 0L, d$y)
cases <- list(
  plain = function() frm(bf(y ~ 1), family = mixture(cumulative(),
                                                      cumulative()),
                         data = d),
  gr = function() frm(bf(y | thres(gr = g) ~ 1),
                      family = mixture(cumulative(), cumulative()), data = d),
  hurdle = function() frm(bf(yh ~ 1),
                          family = mixture(hurdle_cumulative(),
                                           hurdle_cumulative()), data = d),
  hurdle_gr = function() frm(bf(yh2 | thres(gr = g) ~ 1),
                             family = mixture(hurdle_cumulative(),
                                              hurdle_cumulative()),
                             data = d))
for (nm in names(cases)) {
  w <- character(0)
  f <- withCallingHandlers(cases[[nm]](), warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  fam <- f$spec$responses[[1]]$family
  nth <- fam$thres$nthres
  n_id <- sum(nth) + length(nth) * isTRUE(fam$extra_cat)
  cat(sprintf("%-10s nthres %s: parameters %d, lane's identified count %d, Hessian rank %d, warned %s\n",
              nm, paste(nth, collapse = "+"), length(f$opt$par), n_id,
              rank_of(f), any(grepl("no distributional parameter", w))))
}
