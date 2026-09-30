# Reviewer of lane defects, recheck: every reader of the response at the
# rows where an mi() response is missing, on the lane build (fit and
# draws). Log dev/defects-rev-log/mi-frmtmb.txt.
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/defrev-stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
source("dev/defects-rev-mi-data.R")
d <- mi_data()
cap <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning") })
  attr(v, "warnings") <- unique(w)
  v
}
show <- function(k, v) {
  cat(sprintf("  %-18s: ", k))
  if (is.character(v) && length(v) == 1) cat(substr(v, 1, 220)) else
    cat("dim", paste(dim(v) %||% length(v), collapse = "x"), " NA",
        sum(is.na(v)))
  if (length(attr(v, "warnings"))) {
    cat("  [warn:", paste(substr(attr(v, "warnings"), 1, 90), collapse = " | "), "]")
  }
  cat("\n")
}
miss_rows <- list(uni = which(is.na(d$ymi)), mv = which(is.na(d$xmi)),
                  mv2 = which(is.na(d$ymi)), sdy = which(is.na(d$ymeas)))
for (nm in names(mi_models)) {
  cat("=====", nm, " missing rows:", length(miss_rows[[nm]]), "\n")
  f <- cap(frm(eval(mi_models[[nm]]), data = d))
  if (is.character(f)) { cat("  fit:", f, "\n"); next }
  rs <- names(f$spec$responses)
  mr <- if (nm == "mv") "xmi" else rs[1]
  for (ty in c("response", "pearson", "deviance", "osa")) {
    v <- cap(residuals(f, type = ty, resp = mr))
    show(paste("resid", ty), v)
    if (!is.character(v)) {
      e <- if (is.matrix(v)) v[, 1] else v
      cat("     NA exactly at missing rows:",
          identical(which(is.na(e)), miss_rows[[nm]]), "\n")
    }
  }
  show("resid all resp", cap(residuals(f)))
  show("pp_check fit", cap({p <- pp_check(f, resp = mr); "drawn"}))
  show("dharma", cap({dh <- dharma_residuals(f, resp = mr); length(dh$observedResponse)}))
  show("frame y at missing", f$frame$y[[mr]][miss_rows[[nm]]])
  ds <- cap(suppressMessages(frm_sample(f, chains = 1, iter = 400,
                                         refresh = 0, seed = 2)))
  if (is.character(ds)) { cat("  sample:", ds, "\n"); next }
  pe <- cap(predictive_error(ds, resp = mr, ndraws = 20))
  show("predictive_error", pe)
  if (!is.character(pe)) cat("     NA columns exactly the missing rows:",
      identical(which(apply(is.na(pe), 2, all)), miss_rows[[nm]]), "\n")
  show("residuals draws", cap(residuals(ds, resp = mr)))
  r2 <- cap(bayes_R2(ds))
  show("bayes_R2", r2); if (!is.character(r2)) print(r2)
  ll <- cap(log_lik(ds, ndraws = 20))
  show("log_lik", ll)
  if (!is.character(ll)) {
    cat("     log_lik at missing rows (first draw):",
        round(ll[1, miss_rows[[nm]]][1:3], 3), " all finite:",
        all(is.finite(ll)), "\n")
  }
  lo <- cap(loo(ds))
  show("loo", if (is.character(lo)) lo else lo$estimates[1, 1])
  show("pp_check draws", cap({p <- pp_check(ds, resp = mr, ndraws = 5); "drawn"}))
  show("fitted draws lin", cap(fitted(ds, scale = "linear")))
}
