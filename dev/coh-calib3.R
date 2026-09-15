## The false-alarm rate of the shipped assertion.
##
## The assertion is a comparison of two width ratios the run measures:
##
##   ratio = width(cond + (1|id) + (1|id:cond)) / width(cond + (1|id))
##
## in an arm whose truth HAS a subject-by-condition effect and in an arm
## whose truth does not. A test pins one seed per arm, so what has to be
## measured is not "does it pass at that seed" but how far apart the two
## distributions are: an assertion that holds only at a pinned seed has
## broken on this project before.
##
## 16 subjects, 6 rows per cell, sd(id:cond) = 0.5 in the treatment arm.
## dev/coh-calib.R put the two arms at [1.42, 2.60] and [1.00, 1.20]
## over 20 seeds each at this design; this adds 60 more of each.
##
## Run: COH_SEEDS=60 Rscript dev/coh-calib3.R
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
source("dev/coh-sim.R")
source("dev/coh-calib-sim.R")

seeds <- as.integer(Sys.getenv("COH_SEEDS", "60"))
out <- Sys.getenv("COH_OUT", "dev/coh-calib3.tsv")
n_sub <- 16L
n_rep <- 6L
sd_treat <- 0.5

for (arm in c("treat", "null")) {
  sd_ic <- if (identical(arm, "null")) 0 else sd_treat
  for (s in seq_len(seeds)) {
    seed <- 4400000L + s
    d <- coh_flat(seed, n_sub, n_rep, sd_ic)
    t0 <- Sys.time()
    a <- try(coh_fit_one("cond + (1 | id)", d), silent = TRUE)
    b <- try(coh_fit_one("cond + (1 | id) + (1 | id:cond)", d),
             silent = TRUE)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(a, "try-error") || inherits(b, "try-error")) {
      cat("arm=", arm, "\tseed=", seed, "\tERROR\n", sep = "",
          file = out, append = TRUE)
      next
    }
    cat(paste(c(paste0("arm=", arm), paste0("seed=", seed),
                paste0("w_id=", signif(a$hi - a$lo, 8)),
                paste0("w_full=", signif(b$hi - b$lo, 8)),
                paste0("ratio=", signif((b$hi - b$lo) / (a$hi - a$lo),
                                        8)),
                paste0("sd_ic_hat=", signif(b$sd_idcond, 6)),
                paste0("cov_id=", a$lo <= 0.5 && a$hi >= 0.5),
                paste0("cov_full=", b$lo <= 0.5 && b$hi >= 0.5),
                paste0("pd=", a$pdhess && b$pdhess),
                paste0("secs=", signif(secs, 4))),
              collapse = "\t"),
        "\n", sep = "", file = out, append = TRUE)
  }
  message("arm ", arm, " done")
}
