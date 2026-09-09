# lane tmbstan, probe 04: is there a seam that makes the BROKEN build
# constructible in-process? Two candidates, both tried here before
# either is written into the suite.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(testthat)
library(frmtmb)
library(frmtmb.sample)

cat("testthat:", as.character(packageVersion("testthat")), "\n")
cat("local_mocked_bindings exists:",
    "local_mocked_bindings" %in% getNamespaceExports("testthat"), "\n")

det <- frmtmb.sample:::tmbstan_build_broken
env <- environment(det)
cat("closure env has 'cached':", exists("cached", env, inherits = FALSE),
    "\n")

# a model.hpp with the real 2.39 generator shape: TWO log_prob_impl
# overloads, the first patched by autogen and the second left with the
# placeholder. Taken from the shape recorded in
# dev/prior-dropping-investigation.md.
mk_hpp <- function(broken) {
  d <- tempfile("tmbstanhpp"); dir.create(d)
  p <- file.path(d, "model.hpp")
  second <- if (broken) {
    "    lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));"
  } else {
    "    lp_accum__.add(custom_func::custom_func(y));"
  }
  writeLines(c(
    "template <bool propto__, bool jacobian__, typename VecR>",
    "inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__) {",
    "    lp_accum__.add(custom_func::custom_func(y));",
    "    return lp_accum__.sum();",
    "}",
    "template <bool propto__, bool jacobian__, typename VecR>",
    "inline double log_prob_impl(VecR& params_r__) const {",
    second,
    "    return lp_accum__.sum();",
    "}"), p)
  p
}

reset <- function() assign("cached", NULL, envir = env)

cat("\n-- baseline, real installation --\n")
reset(); cat("broken:", det(), "\n")

cat("\n-- candidate A: mock system.file in the package namespace --\n")
res <- tryCatch({
  fake <- mk_hpp(TRUE)
  local({
    testthat::local_mocked_bindings(
      system.file = function(...) fake, .package = "frmtmb.sample")
    reset()
    c(broken = det(),
      errs = inherits(try(frmtmb.sample:::check_tmbstan_build("frm_sample()"),
                          silent = TRUE), "try-error"))
  })
}, error = function(e) paste("FAILED:", conditionMessage(e)))
print(res)

cat("\n-- candidate A, clean arm (the guard must NOT fire) --\n")
res2 <- tryCatch({
  fake <- mk_hpp(FALSE)
  local({
    testthat::local_mocked_bindings(
      system.file = function(...) fake, .package = "frmtmb.sample")
    reset()
    c(broken = det(),
      errs = inherits(try(frmtmb.sample:::check_tmbstan_build("frm_sample()"),
                          silent = TRUE), "try-error"))
  })
}, error = function(e) paste("FAILED:", conditionMessage(e)))
print(res2)

cat("\n-- candidate A, ABSENT arm (no model.hpp at all) --\n")
res3 <- tryCatch({
  local({
    testthat::local_mocked_bindings(
      system.file = function(...) "", .package = "frmtmb.sample")
    reset()
    det()
  })
}, error = function(e) paste("FAILED:", conditionMessage(e)))
print(res3)

cat("\n-- candidate B: poke the cache directly --\n")
assign("cached", TRUE, envir = env)
cat("broken:", det(), "\n")
reset()
cat("after reset, broken:", det(), "\n")

cat("\n-- what class does testthat::skip() signal? --\n")
cnd <- tryCatch(testthat::skip("x"), condition = function(c) c)
cat("classes:", paste(class(cnd), collapse = " / "), "\n")
