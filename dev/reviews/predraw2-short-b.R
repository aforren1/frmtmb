# Short check of predraw2's fix round, item 2: draw_resample()'s k = 0
# for a stratum whose tail weights are all equal, against loo's own k, on
# the equal case and on degenerate ones (one weight left, a cluster of
# equal top weights over a negligible rest).
# Usage (from the lane tree): Rscript dev/reviews/predraw2-short-b.R
local({
.libPaths(c("C:/Users/adf44/source/r/wt-predraw2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), " loo",
    as.character(utils::packageVersion("loo")), "\n")
S <- 1000L
M <- ceiling(min(0.2 * S, 3 * sqrt(S)))
cat("tail size at S = 1000:", M, "\n")
cases <- list(
  all_equal = rep(0, S),
  all_equal_jitter1e16 = rep(0, S) + c(1e-16, -1e-16),
  one_finite_rest_neginf = c(0, rep(-Inf, S - 1)),
  one_top_rest_m1000 = c(0, rep(-1000, S - 1)),
  M_top_rest_m50 = c(rep(0, M), rep(-50, S - M)),
  Mminus1_top_rest_m50 = c(rep(0, M - 1), rep(-50, S - M + 1)),
  half_log03 = rep(c(0, log(0.3)), S / 2),
  all_equal_S24 = rep(0, 24),
  all_equal_S25 = rep(0, 25)
)
for (nm in names(cases)) {
  lw <- cases[[nm]]
  kl <- tryCatch(suppressWarnings(loo::psis(lw, r_eff = 1))$diagnostics$pareto_k,
                 error = function(e) paste("error:", conditionMessage(e)))
  set.seed(3)
  rs <- tryCatch(ns$draw_resample(lw, character(length(lw))),
                 error = function(e) list(k = NA, idx = NA,
                                          err = conditionMessage(e)))
  w <- exp(lw - max(lw)); ess <- sum(w)^2 / sum(w^2)
  cat(sprintf("%-24s n %4d | loo k %s | lane k %s | distinct kept %4d | ESS %.1f %s\n",
              nm, length(lw), format(kl, digits = 3), format(rs$k, digits = 3),
              length(unique(rs$idx)), ess, rs$err %||% ""))
}
# what the lane says on each k (the message or the warning)
cat("\nverdict wording per k (draw_reweight_say, mock fit):\n")
obj <- list(cache = new.env())
for (k in c(0, 0.5, 0.71, Inf)) {
  assign("draw_check_said", FALSE, envir = obj$cache)
  r <- tryCatch(withCallingHandlers(
    ns$draw_reweight_say(obj, list(excess = 2, n = 20L), k, 1000L),
    message = function(m) {
      cat("  k", k, "message:", class(m)[1], "\n"); invokeRestart("muffleMessage")
    }), warning = function(w) cat("  k", k, "warning:", class(w)[1], "\n"))
}
})
