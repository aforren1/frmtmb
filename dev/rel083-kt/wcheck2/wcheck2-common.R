# Lane wcheck2: shared setup for the lane's probe scripts.
#   source("dev/rel083-kt/wcheck2/wcheck2-common.R") after setting `lib_arg` to "base" or a
#   library path. "base" is rellib-r21 (frmtmb 0.82.0).
base_lib <- "C:/Users/adf44/source/r/rellib-r21"
user_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lane_lib <- "C:/Users/adf44/source/r/rellib-r22"
if (!exists("lib_arg")) lib_arg <- commandArgs(TRUE)[1]
.libPaths(if (identical(lib_arg, "base")) c(base_lib, user_lib) else
  c(if (identical(lib_arg, "lane")) lane_lib else lib_arg, user_lib))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")

# One fit's outcome on one line: REFUSED (the weights check), ERROR, or
# the fit's logLik, dispersion coefficients, code and warnings.
outcome <- function(expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, paste(class(w)[1], conditionMessage(w)))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    m <- gsub("\\s+", " ", conditionMessage(r))
    return(paste(if (grepl("has no maximum|is improper", m)) "REFUSED" else
      "ERROR", substr(m, 1, 160)))
  }
  sprintf("fit logLik %.6g betad %s conv %s warn %d%s",
          as.numeric(logLik(r)),
          paste(format(r$estimates$betad, digits = 4), collapse = " "),
          r$opt$convergence, length(ws),
          if (length(ws)) paste0(" [", paste(substr(gsub("\\s+", " ", ws),
                                                    1, 70), collapse = " | "),
                                 "]") else "")
}

# The objective's profile along one dispersion coefficient (index j of
# betad), every other parameter re-optimized from the previous step.
# Returns -fn (the log objective) at each grid value.
profile_disp <- function(ob, g, j = 1L, fix = NULL) {
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  isd <- which(names(p) == "betad")[j]
  fixi <- if (length(fix)) as.integer(names(fix)) else integer()
  if (length(fix)) p[fixi] <- unname(fix)
  free <- setdiff(seq_along(p), c(isd, fixi))
  vapply(g, function(t) {
    p[isd] <- t
    fr <- function(r) {
      q2 <- p
      q2[free] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- tryCatch(nlminb(p[free], fr, control = list(iter.max = 500,
                                                       eval.max = 800)),
                   error = function(e) NULL)
    if (is.null(op)) return(NA_real_)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
}
show_profile <- function(label, g, v) {
  cat(sprintf("%s\n  grid  %s\n  -fn   %s\n  slope %s\n", label,
              paste(g, collapse = " "),
              paste(format(v, digits = 7), collapse = " "),
              paste(format(diff(v) / diff(g), digits = 3), collapse = " ")))
}
q <- function(e) suppressWarnings(suppressMessages(e))
