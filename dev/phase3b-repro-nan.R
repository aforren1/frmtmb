# Reproduce a recovery replicate that stopped at "NA/NaN gradient
# evaluation", and fit the SAME draw two more ways to say whether the
# censoring is what breaks it.
# Usage: Rscript dev/phase3b-repro-nan.R <arm> <seed>
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib2",
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
a <- commandArgs(trailingOnly = TRUE)
arm <- a[[1]]; seed <- as.integer(a[[2]])
# the harness's own simulator, sourced up to its fitting loop
src <- readLines("dev/phase3b-eam-recovery.R")
stop_at <- grep("^for [(]seed in seeds[)]", src)
src <- src[seq_len(stop_at - 1L)]
src <- sub("^args <- commandArgs.*", "args <- c(arm, seed, seed)", src)
eval(parse(text = src))
d <- simulate_one(seed)
try_fit <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(expr, error = function(e) e)
  cat(label, ":", if (inherits(r, "error")) conditionMessage(r) else
    sprintf("ok logLik %.4f code %d", as.numeric(logLik(r)),
            r$opt$convergence),
    sprintf("(%.0f s)", proc.time()[["elapsed"]] - t0), "\n")
  r
}
f1 <- try_fit("as in the harness", fit_one(d))
f_plain <- bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s), ndt ~ 1,
              bias = 0.5)
if (arm == "cens") {
  # the same draw with its response times uncensored: simulate_one()
  # censors only when arm is "cens"
  arm <- "none"
  d0 <- simulate_one(seed)
  f2 <- try_fit("same draw, uncensored, plain wiener",
                frm(f_plain, family = wiener(), data = d0))
} else {
  # the same contaminated draw, plain family
  f2 <- try_fit("same draw, plain wiener, no contaminant",
                frm(f_plain, family = wiener(), data = d))
}
# the failing fit again with core's stage trace
f3 <- try_fit("as in the harness, verbose", {
  arm <- a[[1]]
  if (arm == "cens") {
    frm(bf(rt | dec(upper) + cens(cens) ~ cond + (1 | s),
           bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
        family = wiener(), data = d, verbose = TRUE)
  } else {
    frm(bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s), ndt ~ 1,
           lambda ~ 1, bias = 0.5),
        family = if (arm == "cont") wiener(contaminant = TRUE,
                                           max_ndt = MAXNDT) else
          wiener(contaminant = TRUE),
        data = d, verbose = TRUE)
  }
})
