# Prelude of the latres review scripts. REV_LIB selects the arm:
# "lane" (default) or "base" (rellib-r16 alone).
rev_arm <- Sys.getenv("REV_ARM", "lane")
rev_root <- "C:/Users/adf44/source/r/frmtmb-wt-latres"
user_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
base_lib <- "C:/Users/adf44/source/r/rellib-r16"
lane_lib <- "C:/Users/adf44/source/r/wt-latres-lib"
mut_lib <- "C:/Users/adf44/source/r/rev-latres-lib"
.libPaths(switch(rev_arm,
  lane = c(lane_lib, base_lib, user_lib),
  base = c(base_lib, user_lib),
  mutant = c(mut_lib, lane_lib, base_lib, user_lib)))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(rev_root, "dev/review-latres-out/stan-cache"),
           NOT_CRAN = "true")
Sys.setenv(LANG = "en")
options(warn = 1, width = 110)
setwd(rev_root)
suppressPackageStartupMessages(library(frmtmb))
`%||%` <- function(a, b) if (is.null(a)) b else a
fi <- function(nm) getFromNamespace(nm, "frmtmb")
msg_of <- function(expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(list(val = expr, err = NULL),
             error = function(e) list(val = NULL, err = conditionMessage(e))),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  r$warn <- w
  r
}
