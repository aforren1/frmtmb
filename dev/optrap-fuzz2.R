# Lane optrap, item 2: why last.par.best is not a usable start on fuzz
# seed 20379118. Traces optimizer_from_best() on the lane build: the
# best value the tape recorded, the objective and gradient there now,
# and the run that failed.
#   Rscript dev/optrap-fuzz2.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
env <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-fuzz.R", envir = env)
cat("lib:", dirname(find.package("frmtmb")), "| R:", R.home(), "\n")
sp <- list(family = "Gamma", aterm = "weights", re = "ar1",
           special = "mo_int", dpar = "dpar_x", mode = "reml",
           op = "confint", seed = 20379118L)
d <- env$fuzz_data(sp)
bform <- eval(parse(text = env$fuzz_bf_text(sp)))
ns <- asNamespace("frmtmb")
trace("optimizer_from_best", where = ns, print = FALSE, tracer = quote({
  pb <- obj$env$last.par.best
  p0 <- pb[obj$env$lfixed()]
  cat("  from_best: error:", conditionMessage(e), "\n")
  cat("  value.best recorded:", format(obj$env$value.best, digits = 12),
      "| start of the failed run fn:",
      format(tryCatch(obj$fn(par), error = function(e) NA), digits = 12),
      "\n")
  v <- tryCatch(obj$fn(p0), error = function(e) NA)
  g <- tryCatch(obj$gr(p0), error = function(e) NA)
  cat("  at last.par.best now: fn", format(v, digits = 12),
      "| non-finite gradient entries", sum(!is.finite(g)), "of",
      length(g), "\n")
  cat("  last.par.best outer:", signif(p0, 5), "\n")
}))
r <- tryCatch(suppressWarnings(frm(bform, data = d, REML = TRUE)),
              error = function(e) e)
cat(if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
  sprintf("OK code %d objective %.10f", r$opt$convergence, r$opt$objective),
  "\n")
