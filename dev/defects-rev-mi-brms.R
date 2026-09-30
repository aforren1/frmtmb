# Reviewer of lane defects, recheck: what brms 2.23.0 returns at the rows
# where an mi() response is missing. Log dev/defects-rev-log/mi-brms.txt,
# results dev/defects-rev-log/mi-brms.rds.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
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
out <- list()
for (nm in names(mi_models)) {
  cat("=====", nm, "\n")
  f <- suppressWarnings(brm(eval(mi_models[[nm]]), data = d, chains = 1,
                            iter = 600, refresh = 0, seed = 1))
  resps <- if (is.null(f$formula$responses)) "" else f$formula$responses
  r <- list()
  r$residuals <- cap(residuals(f))
  r$pe <- cap(predictive_error(f, ndraws = 20))
  r$r2 <- cap(bayes_R2(f))
  r$loglik <- cap(log_lik(f, ndraws = 20))
  r$loo <- cap(loo(f))
  r$pp <- cap({p <- pp_check(f, ndraws = 5); ggplot2::ggplot_build(p); "ok"})
  r$fitted_lin <- cap(fitted(f, scale = "linear"))
  for (k in names(r)) {
    v <- r[[k]]
    cat(k, ": ")
    if (is.character(v) && length(v) == 1) cat(substr(v, 1, 200)) else
      cat("dim", paste(dim(v) %||% length(v), collapse = "x"), " NA",
          sum(is.na(v)), " NA cols/rows pattern:",
          if (length(dim(v)) >= 2) paste(which(apply(is.na(v), 1, any))[1:5], collapse = ",") else "")
    if (length(attr(v, "warnings"))) cat("  [warn:", paste(substr(attr(v, "warnings"), 1, 90), collapse = " | "), "]")
    cat("\n")
  }
  if (!is.character(r$r2)) print(r$r2)
  out[[nm]] <- r
}
saveRDS(out, "dev/defects-rev-log/mi-brms.rds")
