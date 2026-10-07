# Lane surface: library paths and environment for a script, by arm.
#
#   source("dev/surface-env.R"); surface_env("lane")   # or "base"
#
# "lane" puts the lane's private library first; "base" measures the
# 0.68.1 reference build alone.
surface_env <- function(arm = c("lane", "base")) {
  arm <- match.arg(arm)
  .libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-surface-lib",
              "C:/Users/adf44/source/r/rellib-r6",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(
    R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
    FRMTMB_STAN_CACHE =
      "C:/Users/adf44/source/r/frmtmb-wt-surface/dev/stan-cache",
    NOT_CRAN = "true")
  invisible(arm)
}

# One line per step: ok or the error, with the warnings counted.
surface_show <- function(tag, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "surface_err"))
  cat(sprintf("%-50s %s%s\n", tag,
              if (inherits(r, "surface_err")) {
                paste("ERROR:", substr(gsub("\n", " ", r), 1, 160))
              } else "ok",
              if (length(w)) {
                paste0(" [", length(w), " warning(s), ", length(unique(w)),
                       " distinct: ", substr(w[1], 1, 70), "]")
              } else ""))
  invisible(r)
}
