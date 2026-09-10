# lane gddm, punch round 1: the review's BLOCKER 2.
#
# gddm_simulate() rep_len()s every parameter to n and then solves once
# per distinct value of `coh`, reading each parameter at that value's
# first trial. A per-trial parameter vector was accepted and half of it
# ignored, with no error and no warning: the same defect the item is
# about, in an exported function.
#
# Seed 31, the review's own, 400 trials, mu = c(-2.5 x 200, +2.5 x 200).
# Arm chosen by GDDM_LIB.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "  frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

n <- 400L
mu <- c(rep(-2.5, n / 2L), rep(2.5, n / 2L))
ctl <- gddm_control(t_max = 2, dt = 0.02, ny = 101L)

run <- function(label, coh) {
  set.seed(31)
  w <- NULL
  r <- withCallingHandlers(
    tryCatch(gddm_simulate(n, mu = mu, bs = 2, ndt = 0.2, coh = coh,
                           control = ctl),
             error = function(e) conditionMessage(e)),
    warning = function(z) { w <<- c(w, conditionMessage(z))
                            invokeRestart("muffleWarning") })
  cat("\n== ", label, "\n", sep = "")
  if (is.character(r)) {
    cat("  refused: ", substr(r, 1, 220), "\n", sep = "")
    return(invisible(NULL))
  }
  h <- seq_len(n / 2L)
  cat(sprintf("  upper rate: first half %.3f  second half %.3f\n",
              mean(r$upper[h]), mean(r$upper[-h])))
  cat("  error or warning:",
      if (is.null(w)) "none" else paste(w, collapse = " / "), "\n")
  invisible(r)
}

run("coh = 0, one condition: mu varies inside it", 0)
run("coh separating the halves", rep(0:1, each = n / 2L))

# The uses that must keep working.
cat("\n== the uses that must keep working\n")
ok <- function(label, expr) {
  set.seed(31)
  r <- tryCatch({
    v <- force(expr)
    sprintf("ok, %d rows, upper rate %.3f", nrow(v), mean(v$upper))
  }, error = function(e) paste("REFUSED:", conditionMessage(e)))
  cat("  ", label, ": ", substr(r, 1, 140), "\n", sep = "")
}
ok("scalar parameters",
   gddm_simulate(100, mu = 1.5, bs = 2, ndt = 0.2, control = ctl))
ok("length-n parameter constant within each coh",
   gddm_simulate(100, mu = rep(c(1, 3), each = 50L), bs = 2, ndt = 0.2,
                 coh = rep(0:1, each = 50L), control = ctl))
ok("a coherence drift, parameters scalar",
   gddm_simulate(100, mu = 2, alpha = 1, bs = 2, ndt = 0.2,
                 coh = rep(c(0, 0.5), each = 50L),
                 drift = gddm_drift_coherence(), control = ctl))
ok("the help page's own example",
   gddm_simulate(20, mu = 1.5, bs = 2, ndt = 0.2,
                 control = gddm_control(t_max = 2, dt = 0.02,
                                        ny = 101)))
