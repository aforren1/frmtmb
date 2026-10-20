# Lane sampleprior, fix round 2: frm_sample() on rdm(), with and
# without ndt_group(), whose run in dev/sampleprior-ndtgroup.R stopped
# with rstan's "no more scalars to read" (data seed 20261121).
#   Rscript dev/sampleprior-rdm.R <lane|r9>
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "r9")) character(0) else
  "C:/Users/adf44/source/r/wt-sampleprior-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(frmtmb.eam)
})
cat("arm", arm, "lib:", find.package("frmtmb.sample"), "\n")
set.seed(20261121)
d <- rdm_simulate(200, v = c(2.4, 1.4), A = 0.5, k = 0.5, ndt = 0.2)
d$subj <- factor(rep(1:4, length.out = nrow(d)))
for (f in list(bf(rt | vint(choice) ~ 1),
               bf(rt | vint(choice) + ndt_group(subj) ~ 1))) {
  fit <- suppressWarnings(frm(f, family = rdm(2), data = d))
  cat("\n==", deparse(f$formula), "\n")
  cat("outer parameters:", length(fit$obj$par), " ML ndt",
      frm_linpred(fit, dpar = "ndt", type = "response")[1L], "\n")
  for (pr in list(NULL, "flat")) {
    r <- tryCatch({
      ds <- suppressWarnings(suppressMessages(frm_sample(
        fit, chains = 2, iter = 300, seed = 4, refresh = 0,
        prior = pr)))
      v <- as.matrix(ds)[, "ndt"]
      sprintf("OK, ndt q01/50/99 %s", paste(sprintf("%.4f",
        stats::quantile(v, c(.01, .5, .99))), collapse = " "))
    }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1,
                                                  60)))
    cat("  prior", if (is.null(pr)) "defaults" else pr, "->", r, "\n")
  }
}
