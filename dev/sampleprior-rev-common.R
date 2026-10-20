# Reviewer of lane sampleprior: shared setup. Captures the objective
# frm_sample() would hand to tmbstan, without sampling, by replacing
# tmbstan::tmbstan in this session only.
lane <- !identical(Sys.getenv("REV_ARM"), "r9")
.libPaths(c(if (lane) "C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")

CAP <- new.env()
real_tmbstan <- tmbstan::tmbstan
capture_on <- function() {
  utils::assignInNamespace("tmbstan", function(obj, ...) {
    CAP$obj <- obj
    CAP$args <- list(...)
    stop("REV_CAPTURED")
  }, "tmbstan")
}
capture_off <- function() {
  utils::assignInNamespace("tmbstan", real_tmbstan, "tmbstan")
}
# the objective frm_sample() samples, from a formula, with `prior`
cap_obj <- function(formula, data, family, prior = NULL, ...) {
  capture_on()
  on.exit(capture_off())
  CAP$obj <- NULL
  r <- tryCatch(suppressMessages(suppressWarnings(
    frm_sample(formula, data = data, family = family, prior = prior,
               chains = 1, iter = 10, refresh = 0, ...))),
    error = function(e) conditionMessage(e))
  if (!identical(r, "REV_CAPTURED") && is.null(CAP$obj)) stop(r)
  CAP$obj
}
# log prior frm_sample() adds at the full parameter vector `p`: the
# default-prior tape less the flat one, both evaluated jointly (env$f
# integrates nothing out)
logprior_at <- function(od, of, p) {
  -(od$env$f(p) - of$env$f(p))
}
