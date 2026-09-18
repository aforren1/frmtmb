source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/stan-cache")
suppressPackageStartupMessages({ library(frmtmb); library(frmtmb.sample) })
cat("tmbstan", format(packageVersion("tmbstan")), "\n")
for (f in c("posterior_epred.frmtmb_draws", "posterior_predict.frmtmb_draws",
            "log_lik.frmtmb_draws", "fitted.frmtmb_fit",
            "predict.frmtmb_fit", "conditional_effects.frmtmb_fit")) {
  g <- tryCatch(get(f, envir = asNamespace(
        if (grepl("draws", f)) "frmtmb.sample" else "frmtmb")),
        error = function(e) NULL)
  cat(sprintf("%-34s allow_new_levels in formals: %s\n", f,
              if (is.null(g)) "NOT FOUND" else
                "allow_new_levels" %in% names(formals(g))))
}
cat("\nbrms's own:\n")
suppressPackageStartupMessages(library(brms))
for (f in c("posterior_epred.brmsfit", "posterior_predict.brmsfit",
            "log_lik.brmsfit", "fitted.brmsfit", "predict.brmsfit")) {
  g <- get(f, envir = asNamespace("brms"))
  cat(sprintf("%-28s allow_new_levels: %s\n", f,
              "allow_new_levels" %in% names(formals(g))))
}
