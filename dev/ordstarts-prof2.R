# Lane ordstarts, build fix round 1: what the objective test of
# mixture_swap_distinct() costs on a mixture of one family with random
# effects (the swap is the same model there, so one evaluation is all it
# adds). Per seed, in one process: the fit with starts = 0 and the
# default, interleaved, and the default's starts record.
#   Rscript dev/ordstarts-prof2.R <family> <re> <n> <G> <seeds>
a <- commandArgs(trailingOnly = TRUE)
famname <- a[1]
re <- a[2]
nobs <- as.integer(a[3])
G <- as.integer(a[4])
seeds <- eval(parse(text = a[5]))
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
run <- readLines("dev/ordstarts-run.R")
eval(parse(text = run[grep("^mkg <- ", run):(grep("^form <- ", run) - 1L)]))
form <- switch(re, int = bf(y ~ x + z + (1 | g)),
               slope = bf(y ~ x + z + (1 + x | g)))
fam <- switch(famname, cum_cum = mixture(cumulative(), cumulative()),
              cum_sratio = mixture(cumulative(), sratio()))
q <- function(e) suppressWarnings(suppressMessages(e))
r <- NULL
for (s in seeds) {
  d <- mkg(s)
  t0 <- proc.time()[["elapsed"]]
  f0 <- tryCatch(q(frm(form, family = fam, data = d,
                       control = frmtmb_control(starts = 0))),
                 error = function(e) NULL)
  e0 <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  f1 <- tryCatch(q(frm(form, family = fam, data = d)),
                 error = function(e) NULL)
  e1 <- proc.time()[["elapsed"]] - t0
  if (is.null(f0) || is.null(f1)) next
  ms <- f1$opt$mixture_starts
  r <- rbind(r, data.frame(seed = s, e0 = e0, e1 = e1, ratio = e1 / e0,
                           same_ll = abs(as.numeric(logLik(f1)) -
                                           as.numeric(logLik(f0))) < 1e-9,
                           rec = if (is.null(ms)) "NULL" else
                             paste(ms, collapse = ",")))
}
print(r, digits = 3)
cat("PROF2", famname, re, nobs, "ratio median", signif(median(r$ratio), 3),
    "sum", signif(sum(r$e1) / sum(r$e0), 3), "\n")
