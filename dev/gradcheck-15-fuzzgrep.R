# Punch item 2 kept the warning's opening phrase because
# tests/testthat/helper-fuzz.R greps it in FUZZ_NONCONVERGENCE. A gated
# fuzz run that passes does not PROVE the grep still matches, because it
# passes either way when nothing fires. This asserts the match directly,
# on every branch of the rewritten message, with the pattern read out of
# the harness source rather than retyped.
#
#   Rscript dev/gradcheck-15-fuzzgrep.R lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "lane"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

# the harness's own pattern, taken from the harness
env <- new.env()
sys.source("tests/testthat/helper-fuzz.R", envir = env)
pat <- get("FUZZ_NONCONVERGENCE", envir = env)
cat("FUZZ_NONCONVERGENCE from helper-fuzz.R:\n  ", pat, "\n\n")

cases <- list(
  bounded_short = list(gmax = 5971.25, gmax_par = "x", proj = 52.4887,
                       proj_par = "z", bound_held = "x",
                       headroom = 0.130833, tol = 1e-3),
  plain_short = list(gmax = 0.004986, gmax_par = "(Intercept)",
                     proj = 0.004986, proj_par = "(Intercept)",
                     bound_held = character(0), headroom = 0.5295,
                     tol = 1e-3),
  unusable = list(gmax = 50170, gmax_par = "sigma_(Intercept)",
                  proj = 50170, proj_par = "sigma_(Intercept)",
                  bound_held = character(0), headroom = NA_real_,
                  tol = 1e-3),
  two_bounds = list(gmax = 147.1, gmax_par = "x1", proj = 2.58e-07,
                    proj_par = "(Intercept)",
                    bound_held = c("x1", "x2"), headroom = 1e-09,
                    tol = 1e-3),
  no_names = list(gmax = 12.5, gmax_par = NA_character_, proj = 12.5,
                  proj_par = NA_character_, bound_held = character(0),
                  headroom = 0.02, tol = 1e-3))

ok <- TRUE
for (nm in names(cases)) {
  m <- frmtmb:::grad_warning_msg(cases[[nm]])
  hit <- grepl(pat, m)
  ok <- ok && hit
  cat(sprintf("%-14s grep %-5s  %s\n", nm, hit, m))
  cat("\n")
}
cat("every branch matched:", ok, "\n")
if (!ok) stop("the fuzz harness would stop recognizing this warning")
