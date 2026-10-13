# Reviewer of the 0.70.0 consolidation: on the cox() study's data
# (dev/moparity-cox-study.R, n = 300, seeds 1..N), the fits whose
# mo_search() refits nothing (runs 0) but probes a face, on the release
# build and on lane moparity's: the objective's state after the fit
# (last and best point at the fitted optimum) and the standard errors,
# seed by seed. The merge dropped moparity's restore after the probes in
# favor of optrap's obj_settle(); this shows what the fit reads after.
#   Rscript dev/relrev070-cox-runs.R r8|lane N
args <- commandArgs(TRUE)
arm <- args[1]
N <- as.integer(args[2])
.libPaths(c(switch(arm, r8 = "C:/Users/adf44/source/r/rellib-r8",
                   lane = "C:/Users/adf44/source/r/wt-moparity-lib"),
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, as.character(packageVersion("frmtmb")),
    find.package("frmtmb"), "\n")
state_at <- function(obj, par) {
  env <- obj$env
  lf <- env$lfixed()
  at <- function(p) identical(unname(as.numeric(p[lf])),
                              unname(as.numeric(par)))
  at(env$last.par) && at(env$last.par.best)
}
rows <- list()
for (seed in seq_len(N)) {
  set.seed(seed)
  n <- 300
  x <- rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
  fit <- suppressWarnings(frm(bf(time | cens(cens) ~ x), family = cox(),
                              data = dd))
  s <- as.numeric(cox_baseline(fit))
  ms <- fit$opt$mo_search
  at <- state_at(fit$obj, fit$opt$par)
  se <- summary(fit)$fixed[, 2]
  rows[[seed]] <- data.frame(
    seed = seed, runs = if (is.null(ms)) 0 else ms[["runs"]],
    polish = if (is.null(fit$opt$polish)) 0 else fit$opt$polish[["steps"]],
    minw = min(s), state_at = at, ll = as.numeric(logLik(fit)),
    se_x = unname(se[2]), se_int = unname(se[1]),
    nse_nan = sum(!is.finite(se)))
}
res <- do.call(rbind, rows)
saveRDS(res, file.path("dev/relrev070-log", paste0("cox-runs-", arm, ".rds")))
probe0 <- res$runs == 0 & res$minw < 1e-6
cat("seeds", nrow(res), "| runs 0:", sum(res$runs == 0),
    "| runs 0 with a face (probed):", sum(probe0),
    "| polished:", sum(res$polish > 0),
    "| state at the optimum after the fit:", sum(res$state_at),
    "of", nrow(res), "| probed and at:", sum(res$state_at[probe0]),
    "| any NaN SE:", sum(res$nse_nan > 0), "\n")
