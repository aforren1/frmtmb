## Probe: does dev/coh-sim.R draw the tier's data, and what do the
## accessors call the condition contrast?
##
## Run:
##   Rscript dev/coh-probe.R
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
source("dev/coh-sim.R")

cat("frmtmb", format(packageVersion("frmtmb")), " coupling",
    format(packageVersion("frmtmb.coupling")), "\n")

## 1. the data, against the tier's own generator, byte for byte
tier <- new.env()
tier$scale_small <- function() FALSE
src <- readLines("extensions/frmtmb.coupling/tests/testthat/test-scale.R")
stop_at <- grep("^for [(]m in coupling_models", src)[1L]
eval(parse(text = paste(src[seq_len(stop_at - 1L)], collapse = "\n")),
     envir = tier)
a <- tier$coupling_scale_data(20260908L)
b <- coh_data(20260908L)
cat("data identical to the tier:", identical(a, b), "\n")

## 2. one fit of the correct rung, against dev/scale-findings.md
r <- coh_fit_one(coh_rungs()$full, b)
cat(sprintf(paste0("full: est=%.4f se=%.4f (%.3f, %.3f) sd_id=%.4f ",
                   "sd_idcond=%.4f ll=%.3f conv=%s pd=%s %.1fs\n"),
            r$est, r$se, r$lo, r$hi, r$sd_id, r$sd_idcond, r$loglik,
            r$conv, r$pdhess, r$secs))
r2 <- coh_fit_one(coh_rungs()$id, b)
cat(sprintf(paste0("id:   est=%.4f se=%.4f (%.3f, %.3f) sd_id=%.4f ",
                   "ll=%.3f conv=%s pd=%s %.1fs\n"),
            r2$est, r2$se, r2$lo, r2$hi, r2$sd_id, r2$loglik,
            r2$conv, r2$pdhess, r2$secs))
cat("width ratio id / full:", (r2$hi - r2$lo) / (r$hi - r$lo), "\n")
