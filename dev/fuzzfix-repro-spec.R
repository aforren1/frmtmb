# Rebuild one of the two fuzz specs from the generator and run it.
#   Rscript dev/fuzzfix-repro-spec.R <lib|base> <cum|gam>
# <lib> "base" means rellib-r10 alone.
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
source("tests/testthat/helper-fuzz.R")
sp <- switch(a[2],
  cum = list(family = "cumulative", aterm = "weights", re = "nested",
             special = "none", dpar = "none", mode = "quadrature",
             op = "vcov", seed = 20322452L, kind = "cover"),
  gam = list(family = "Gamma", aterm = "weights", re = "ar1",
             special = "mo_int", dpar = "dpar_x", mode = "profile",
             op = "vcov", seed = 20502220L, kind = "pad"))
d <- fuzz_data(sp)
cat("n:", nrow(d), "\n")
cat(fuzz_call_text(sp), "\n")
saveRDS(d, sprintf("dev/fuzzfix-out/data-%s.rds", a[2]))
r <- fuzz_fit_one(sp, d)
cat("ok:", r$ok, "\n")
if (!r$ok) print(conditionMessage(r$error))
for (w in r$warnings) cat("WARN:", w, "\n")
if (r$ok) {
  fit <- r$value
  print(logLik(fit))
  s <- tryCatch(withCallingHandlers(summary(fit),
    error = function(e) print(sys.calls())), error = function(e) e)
  if (inherits(s, "error")) cat("SUMMARY ERROR:", conditionMessage(s), "\n")
  else print(s)
}
cat("END\n")
