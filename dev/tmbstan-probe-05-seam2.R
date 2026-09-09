# lane tmbstan, probe 05: two more seams for the BROKEN construction.
# C: shadow system.file in the closure's own local() environment.
# D: prepend a synthesized tmbstan to .libPaths() while tmbstan is not
#    loaded, which is the seam that needs no knowledge of the closure.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
library(frmtmb.sample)

det <- frmtmb.sample:::tmbstan_build_broken
env <- environment(det)
cat("closure env is the namespace:", identical(env,
    asNamespace("frmtmb.sample")), "\n")
cat("closure env locked:", environmentIsLocked(env), "\n")

mk_hpp <- function(broken) {
  d <- tempfile("hpp"); dir.create(d)
  p <- file.path(d, "model.hpp")
  second <- if (broken) {
    "    lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));"
  } else {
    "    lp_accum__.add(custom_func::custom_func(y));"
  }
  writeLines(c("inline auto log_prob_impl(VecR& p) {",
               "    lp_accum__.add(custom_func::custom_func(y));",
               "}",
               "inline double log_prob_impl(VecR& p) const {",
               second, "}"), p)
  p
}
reset <- function() assign("cached", NULL, envir = env)

cat("\n== seam C: shadow system.file inside the closure env ==\n")
for (br in c(TRUE, FALSE)) {
  p <- mk_hpp(br)
  assign("system.file", function(...) p, envir = env)
  reset()
  d <- det()
  e <- tryCatch({
    frmtmb.sample:::check_tmbstan_build("frm_sample()"); NA_character_
  }, error = function(e) conditionMessage(e))
  cat("broken hpp =", br, "-> detector =", d,
      "| refusal =", if (is.na(e)) "silent" else "THROWS", "\n")
  if (!is.na(e)) cat("   message starts: ", substr(e, 1, 72), "\n")
}
# the ABSENT arm: no model.hpp at all
assign("system.file", function(...) "", envir = env)
reset()
cat("no model.hpp at all -> detector =", det(), "\n")
rm("system.file", envir = env)
reset()
cat("after removing the shadow, detector =", det(), " (real install)\n")

cat("\n== seam D: a synthesized tmbstan on .libPaths() ==\n")
cat("tmbstan loaded at this point:",
    "tmbstan" %in% loadedNamespaces(), "\n")
tl <- tempfile("lanelib"); dir.create(file.path(tl, "tmbstan"),
                                      recursive = TRUE)
writeLines(c("Package: tmbstan", "Version: 1.2.0",
             "Built: R 4.6.1; ; 2026-01-01 00:00:00 UTC; windows"),
           file.path(tl, "tmbstan", "DESCRIPTION"))
file.copy(mk_hpp(TRUE), file.path(tl, "tmbstan", "model.hpp"))
old <- .libPaths()
.libPaths(c(tl, old))
cat("system.file resolves to:",
    system.file("model.hpp", package = "tmbstan"), "\n")
reset()
cat("detector =", det(), "\n")
.libPaths(old)
reset()
cat("restored, detector =", det(), "\n")

cat("\n== seam D under a LOADED tmbstan, which is the suite's state ==\n")
loadNamespace("tmbstan")
.libPaths(c(tl, old))
cat("system.file resolves to:",
    system.file("model.hpp", package = "tmbstan"), "\n")
reset()
cat("detector =", det(), "\n")
.libPaths(old)
