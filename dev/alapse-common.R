# Lane alapse: shared setup for the lane's probe scripts.
#   source("dev/alapse-common.R") after setting `lib_arg` to "base"
#   (rellib-r22, frmtmb 0.83.0) or "lane" (wt-alapse-lib before it).
base_lib <- "C:/Users/adf44/source/r/rellib-r22"
user_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lane_lib <- "C:/Users/adf44/source/r/wt-alapse-lib"
if (!exists("lib_arg")) lib_arg <- commandArgs(TRUE)[1]
.libPaths(if (identical(lib_arg, "base")) c(base_lib, user_lib) else
  c(lane_lib, base_lib, user_lib))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")

# Run expr, keeping every warning and message with its class.
cap <- function(expr) {
  w <- character(); wc <- character(); m <- character(); mc <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(z) {
      w <<- c(w, gsub("\\s+", " ", conditionMessage(z)))
      wc <<- c(wc, class(z)[1])
      invokeRestart("muffleWarning")
    },
    message = function(z) {
      m <<- c(m, gsub("\\s+", " ", conditionMessage(z)))
      mc <<- c(mc, class(z)[1])
      invokeRestart("muffleMessage")
    })
  list(v = v, w = w, wc = wc, m = m, mc = mc)
}

# The asymmetric Laplace log likelihood, with weights.
ald_ll <- function(y, mu, s, p, w = 1) {
  r <- y - mu
  sum(w * (log(p) + log1p(-p) - log(s) - r * (p - (r < 0)) / s))
}
