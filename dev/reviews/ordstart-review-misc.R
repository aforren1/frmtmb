# Review of lane ordstart: odds and ends, per arm (base | lane).
# 1. frm_allfit(start_from_mle = FALSE) time on cum_sratio seed 4;
# 2. profile confint() on one parameter (does it run the starts?);
# 3. cum_cum seed 10 with x rescaled by 1e-4 (autoscale): the degenerate
#    warning with and without the rescaling.
#   Rscript ordstart-review-misc.R arm
arm <- commandArgs(trailingOnly = TRUE)[1]
libs <- switch(arm, base = "C:/Users/adf44/source/r/rellib-r20",
               lane = c("C:/Users/adf44/source/r/wt-ordstart-lib",
                        "C:/Users/adf44/source/r/rellib-r20"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
q <- function(e) suppressWarnings(suppressMessages(e))
tm <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(q(eval.parent(substitute(expr))), error = function(e) e)
  cat(sprintf("TIME %-36s %7.2f s %s\n", label,
              proc.time()[["elapsed"]] - t0,
              if (inherits(r, "error")) conditionMessage(r) else ""))
  invisible(r)
}
d <- mk(4, mix = TRUE)
f <- q(frm(bf(y ~ x + z), family = mixture(cumulative(), sratio()), data = d))
tm("frm_allfit(start_from_mle = FALSE)", af <- frm_allfit(f, start_from_mle = FALSE))
print(af)
tm("confint profile mu1_x", confint(f, parm = "mu1_x", method = "profile"))
d <- mk(10, mix = TRUE)
d$xs <- d$x * 1e-4
for (fo in list(bf(y ~ x + z), bf(y ~ xs + z))) {
  w <- character()
  g <- withCallingHandlers(frm(fo, family = mixture(cumulative(), cumulative()),
                               data = d),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  cat(sprintf("SEED10 %s ll %.4f autoscaled %s degen %s\n",
              deparse(fo$formula), as.numeric(logLik(g)), isTRUE(g$autoscaled),
              any(grepl("ended at a degenerate boundary", w))))
}
cat("DONE\n")
