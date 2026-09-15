# Lane `learnhier`, Test 1: is the `sd(ndt)` shrinkage worse for
# learners whose conditional is more sharply bounded?
#
#   Rscript dev/learnhier-test1.R
#
# NO NEW FITS. It reads the per-learner vectors that
# `dev/learnhier-run2.R` saves and that the original runner did not.
#
# THE PREDICTION WAS WRITTEN DOWN FIRST, in dev/learnhier-findings.md
# under "Test 1's prediction". The Laplace error at the hard edge of the
# Wiener density says the attenuation of a learner's fitted deviation is
# STRONGEST where its observed floor sits CLOSEST to its own true
# non-decision time, so the slope of the fitted deviation on the true
# one RISES with `margin`. The competing explanations, the ordinary
# downward bias of a maximum likelihood variance component and boundary
# attraction on a correlated pair, are properties of the BLOCK and
# predict no relationship with an individual learner's margin at all. A
# flat slope therefore refutes the last candidate standing.
#
# WITHIN REPLICATE, NOT BETWEEN. Replicates differ in their realized
# spreads and in their floors, so `u_true`, `u_hat` and `margin` are
# each standardized inside their own replicate before anything is
# pooled. A between-replicate association would not be evidence for a
# claim about individual learners.
#
# THE CONTROL. `drift` has no bound on its conditional at all, so the
# same regression run on the `drift` component of the same block must
# come back flat. If it does not, the statistic is picking up something
# that is not the edge, and the `ndt` result cannot be read.
source("dev/learnhier-env.R")

fs <- list.files("dev/learnhier-rec", pattern = "^rlddm-[0-9]+[.]rds$",
                 full.names = TRUE)
R <- Filter(function(r) !is.null(r$per_learner),
            lapply(fs, readRDS))
cat("replicates carrying per-learner vectors: ", length(R), " of ",
    length(fs), "\n", sep = "")
if (!length(R)) {
  cat("none yet; dev/learnhier-run2.R writes them\n")
  quit(save = "no")
}

zs <- function(x) (x - mean(x)) / stats::sd(x)

# One long frame: every learner of every replicate, standardized inside
# its own replicate.
mk <- function(col) {
  do.call(rbind, lapply(seq_along(R), function(i) {
    p <- R[[i]]$per_learner
    ut <- p$dev_true[, col]
    uh <- p$ranef[, col]
    data.frame(rep = i, seed = R[[i]]$seed,
               u_true = zs(ut), u_hat = (uh - mean(uh)) / stats::sd(ut),
               margin = zs(p$margin), margin_raw = p$margin,
               stringsAsFactors = FALSE)
  }))
}

report <- function(d, label, control) {
  cat("\n== ", label, if (control) "  (CONTROL, no bound)" else "",
      " ==\n", sep = "")
  cat("learners ", nrow(d), " over ", length(unique(d$rep)),
      " replicates\n", sep = "")
  # the overall attenuation, for context: a slope of 1 is no shrinkage
  m0 <- stats::lm(u_hat ~ 0 + u_true, data = d)
  cat("overall slope of the fitted deviation on the true one: ",
      round(stats::coef(m0)[["u_true"]], 4), "\n", sep = "")
  # the test: does that slope depend on the margin?
  m1 <- stats::lm(u_hat ~ 0 + u_true + u_true:margin, data = d)
  ci <- stats::confint(m1)
  k <- "u_true:margin"
  cat("interaction slope x margin: ",
      round(stats::coef(m1)[[k]], 4), "  95% CI (",
      round(ci[k, 1L], 4), ", ", round(ci[k, 2L], 4), ")  p = ",
      format(summary(m1)$coefficients[k, 4L], digits = 3), "\n", sep = "")
  # and the same thing without a model, by margin quartile within
  # replicate, because a reader should not have to trust one regression
  d$q <- unlist(lapply(split(d$margin_raw, d$rep), function(x) {
    as.integer(cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE))
  }))
  qs <- vapply(1:4, function(q) {
    s <- d[d$q == q, , drop = FALSE]
    stats::coef(stats::lm(u_hat ~ 0 + u_true, data = s))[["u_true"]]
  }, numeric(1))
  cat("slope by margin quartile within replicate, narrowest first: ",
      paste(round(qs, 4), collapse = "  "), "\n", sep = "")
  cat("  quartile 4 minus quartile 1: ", round(qs[[4L]] - qs[[1L]], 4),
      "\n", sep = "")
  invisible(list(inter = stats::coef(m1)[[k]], ci = ci[k, ], q = qs,
                 slope = stats::coef(m0)[["u_true"]]))
}

cols <- colnames(R[[1L]]$per_learner$dev_true)
res <- list()
for (k in cols) res[[k]] <- report(mk(match(k, cols)), k, k != "ndt")

# THE CONTROLS DID NOT COME BACK FLAT, and that is a property of the
# STATISTIC rather than of the model. `margin_i = floor_i - ndt_i` is
# the smallest first-passage time that learner produced, so it is a
# function of that learner's OWN drift and boundary: a faster learner
# has a smaller minimum passage time and therefore a smaller margin. The
# margin is correlated with the other deviations in the block, and any
# shrinkage that varies with a learner's speed shows up as a margin
# interaction for EVERY component, bounded conditional or not.
#
# So the raw ndt interaction cannot be read on its own. What can be read
# is the EXCESS over the unbounded components, taken per replicate so
# both come from the same learners and their errors are paired. This
# adjustment is made because the CONTROL FAILED, which is an
# instrument-driven reason, and not because of anything the ndt number
# said. The raw figures are printed above either way.
cat("
== the controls did not come back flat; the paired excess ==
")
per_rep <- function(col) {
  d <- mk(match(col, cols))
  vapply(split(d, d$rep), function(s) {
    stats::coef(stats::lm(u_hat ~ 0 + u_true + u_true:margin,
                          data = s))[["u_true:margin"]]
  }, numeric(1))
}
pr <- vapply(cols, per_rep, numeric(length(R)))
ctrl <- rowMeans(pr[, cols != "ndt", drop = FALSE])
dif <- pr[, "ndt"] - ctrl
tt <- stats::t.test(dif)
cat("ndt interaction minus the mean of the three unbounded ones,",
    " per replicate
", sep = "")
cat("  n = ", length(dif), " replicates, mean ", round(mean(dif), 4),
    ", 95% CI (", round(tt$conf.int[[1L]], 4), ", ",
    round(tt$conf.int[[2L]], 4), "), p = ",
    format(tt$p.value, digits = 3), "
", sep = "")
cat("  positive on ", sum(dif > 0), " of ", length(dif),
    " replicates
", sep = "")
cat("
attenuation by component (a slope of 1 is no shrinkage):
")
for (k in cols) {
  cat("  ", formatC(k, width = 6), " ", round(res[[k]]$slope, 4), "
",
      sep = "")
}
cat("The gap any margin story has to explain is ",
    round(res[["drift"]]$slope - res[["ndt"]]$slope, 3),
    ", ndt against drift.
", sep = "")
