# Reviewer probe for lane tmbstan: construct the case where the guarded
# thing is ABSENT rather than merely wrong.
#
# Three absences, through the same seam the lane's own fixture uses
# (a `system.file` bound in the detector's local() environment):
#   1. system.file() answers "" , the shape of a tmbstan that ships no
#      model.hpp;
#   2. model.hpp exists but is empty;
#   3. model.hpp carries the defect under a marker string the detector
#      does not know, which is what an upstream rename would look like.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
library(testthat)

env <- environment(frmtmb.sample:::tmbstan_build_broken)

with_hpp <- function(answer, lines) {
  if (!is.null(lines)) writeLines(lines, answer)
  old_exists <- exists("system.file", envir = env, inherits = FALSE)
  assign("system.file", function(...) answer, envir = env)
  assign("cached", NULL, envir = env)
  on.exit({
    if (!old_exists) rm("system.file", envir = env)
    assign("cached", NULL, envir = env)
  }, add = TRUE)
  broken <- frmtmb.sample:::tmbstan_build_broken()
  refused <- inherits(try(frmtmb.sample:::check_tmbstan_build("frm_sample()"),
                          silent = TRUE), "try-error")
  skipped <- tryCatch({
    testthat::skip_if_not_installed("tmbstan")
    if (frmtmb.sample:::tmbstan_build_broken()) testthat::skip("x")
    "ran"
  }, skip = function(cnd) "skipped")
  c(broken = broken, refused = refused, gate = skipped)
}

tf <- tempfile(fileext = ".hpp")
row <- function(label, ...) {
  r <- with_hpp(...)
  cat(sprintf("%-34s detector=%-5s refuses=%-5s skip_sampler=%s\n",
              label, r[["broken"]], r[["refused"]], r[["gate"]]))
}

row("1 no model.hpp at all", "", NULL)
row("2 model.hpp present but empty", tf, character(0))
row("3 defect under a renamed marker", tf, c(
  "inline auto log_prob_impl(VecR& params_r__) {",
  "    lp_accum__.add(custom_func::custom_func(y));",
  "}",
  "inline double log_prob_impl(VecR& params_r__) const {",
  "    lp_accum__.add(stan::math::std_normal_lpdf<propto__, false>(y));",
  "}"))
row("4 the defect as it ships today", tf, c(
  "inline auto log_prob_impl(VecR& params_r__) {",
  "    lp_accum__.add(custom_func::custom_func(y));",
  "}",
  "inline double log_prob_impl(VecR& params_r__) const {",
  "    lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));",
  "}"))

cat("\nstate handed back: detector =",
    frmtmb.sample:::tmbstan_build_broken(), "\n")
