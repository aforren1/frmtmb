# Lane optrap, item 2: fuzz seed 20379118 (Gamma, REML, weights,
# ar1(tt + 0 | g), mo(xo) * z, shape ~ 1 + z), the spec
# dev/relrev069-fuzzdiag.R traces. Fits it and reports the outcome, and
# where the objective's best point is.
#   Rscript dev/optrap-fuzz1.R <lib or ""> [seed]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seed <- if (length(a) > 1) as.integer(a[2]) else 20379118L
suppressMessages({library(testthat); library(frmtmb)})
env <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-fuzz.R", envir = env)
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "| R:", R.home(), "\n")
sp <- list(family = "Gamma", aterm = "weights", re = "ar1",
           special = "mo_int", dpar = "dpar_x", mode = "reml",
           op = "confint", seed = seed)
d <- env$fuzz_data(sp)
bform <- eval(parse(text = env$fuzz_bf_text(sp)))
w <- character()
r <- tryCatch(withCallingHandlers(
  frm(bform, data = d, REML = TRUE, verbose = TRUE),
  warning = function(x) {w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")},
  message = function(x) {cat("  ", conditionMessage(x))
    invokeRestart("muffleMessage")}), error = function(e) e)
if (inherits(r, "error")) {
  cat("ERROR:", substr(conditionMessage(r), 1, 200), "\n")
} else {
  e <- r$obj$env
  lf <- e$lfixed()
  cat(sprintf("OK code %d (%s) objective %.10f\n", r$opt$convergence,
              r$opt$message, r$opt$objective))
  cat("best point is the optimum:",
      identical(unname(e$last.par.best[lf]), unname(r$opt$par)),
      "| last point is the optimum:",
      identical(unname(e$last.par[lf]), unname(r$opt$par)), "\n")
  cat("finite gradient at the optimum:",
      all(is.finite(r$obj$gr(r$opt$par))), "\n")
  for (x in w) cat("  W:", substr(x, 1, 150), "\n")
}
