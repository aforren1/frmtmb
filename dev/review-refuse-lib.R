# Shared setup of the refuse review's scripts. Run from the worktree root.
.libPaths(c(strsplit(Sys.getenv("REV_LIB", "C:/Users/adf44/source/r/wt-refuse-lib"),
                     ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-refuse-out/stan-cache",
                                             winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
if (!"equiv_brms_me" %in% search()) {
  attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
}
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
# one spec through the survey's own flow
meas <- function(spec) {
  E$eq_specs[["variant"]] <- function() spec
  rec <- E$eq_measure("variant", list(status = "untested"))
  keep <- c("verdict", "reason", "frm_error", "brms_error", "l1_const",
            "l1_grad", "l1_extra_const", "l1_error")
  rec[intersect(names(rec), keep)]
}
pr <- function(...) cat(sprintf(...), "\n", sep = "")
tr <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(list(msg = conditionMessage(e)), class = "tr_err")
    }),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  attr(v, "w") <- unique(w)
  v
}
is_err <- function(v) inherits(v, "tr_err")
msg <- function(v, n = 200) {
  if (is_err(v)) substr(gsub("[\r\n]+", " ", v$msg), 1, n) else "ok"
}
