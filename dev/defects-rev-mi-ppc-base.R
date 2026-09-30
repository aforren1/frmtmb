# Reviewer of lane defects, recheck: dev/defects-rev-mi-ppc.R on the base build.
# response, on a draws object of an mi() fit (lane build).
.libPaths(c(
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/defrev-stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
source("dev/defects-rev-mi-data.R")
d <- mi_data()
f <- frm(ymi | mi() ~ x, data = d)
ds <- suppressWarnings(suppressMessages(frm_sample(f, chains = 1, iter = 400, refresh = 0, seed = 2)))
for (ty in c("dens_overlay", "error_hist", "error_scatter_avg", "scatter_avg", "stat", "intervals", "loo_pit_overlay")) {
  w <- character(0)
  r <- withCallingHandlers(tryCatch({p <- pp_check(ds, type = ty, ndraws = 5); ggplot2::ggplot_build(p); "ok"},
         error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 150))),
       warning = function(cnd) { w <<- c(w, substr(conditionMessage(cnd), 1, 80)); invokeRestart("muffleWarning") })
  cat(sprintf("%-18s %s  [warn: %s]\n", ty, r, paste(unique(w), collapse = " | ")))
}
if (exists("loo_R2")) print(tryCatch(loo_R2(ds), error = function(e) conditionMessage(e)))
