## The shipped assertion, calibrated against a reference the RUN
## measures rather than against a constant.
##
## dev/coh-calib.R showed that comparing the two arms' width ratios at
## one seed each does not separate them: at 16 subjects the treatment
## arm reaches down to 1.09 and the null arm up to 1.38. So the test
## needs a reference for what the interval SHOULD be, computed from the
## same data and without a model.
##
## The two-stage reference. Pool each subject-by-condition cell into one
## Hermitian matrix, take its magnitude squared coherence, put it on the
## logit scale, difference the two conditions within a subject and treat
## the
## per-subject differences as a sample:
##
##   se_ref = sd(d_i) / sqrt(n_sub)
##
## It is inefficient and its point estimate is attenuated by the
## upward bias of a coherence from finitely many segments, which is why
## nothing here compares ESTIMATES. What it is, is a standard error for
## a within-subject contrast that carries the subject-by-condition
## spread by construction, because that spread is in the d_i.
##
## Run: COH_SEEDS=40 Rscript dev/coh-calib2.R
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
source("dev/coh-sim.R")
source("dev/coh-calib-sim.R")

seeds <- as.integer(Sys.getenv("COH_SEEDS", "40"))
out <- Sys.getenv("COH_OUT", "dev/coh-calib2.tsv")
## Three designs rather than dev/coh-calib.R's eight: that sweep
## already showed the width ratio moves the same way with every one of
## the three knobs, so what is left is whether the reference separates
## the arms at the CHEAPEST of them.
grid <- data.frame(n_sub = c(16L, 16L, 24L), n_rep = c(6L, 6L, 10L),
                   sd_ic = c(0.35, 0.5, 0.35))

for (g in seq_len(nrow(grid))) {
  for (arm in c("treat", "null")) {
    sd_ic <- if (identical(arm, "null")) 0 else grid$sd_ic[g]
    for (s in seq_len(seeds)) {
      d <- coh_flat(4300000L + 1000L * g + s, grid$n_sub[g],
                    grid$n_rep[g], sd_ic)
      ref <- coh_two_stage(d)
      a <- try(coh_fit_one("cond + (1 | id)", d), silent = TRUE)
      b <- try(coh_fit_one("cond + (1 | id) + (1 | id:cond)", d),
               silent = TRUE)
      if (inherits(a, "try-error") || inherits(b, "try-error")) next
      cat(paste(c(paste0("g=", g), paste0("n_sub=", grid$n_sub[g]),
                  paste0("n_rep=", grid$n_rep[g]),
                  paste0("sd_ic_true=", sd_ic), paste0("arm=", arm),
                  paste0("s=", s),
                  paste0("se_id=", signif(a$se, 8)),
                  paste0("se_full=", signif(b$se, 8)),
                  paste0("se_ref=", signif(ref$se, 8)),
                  paste0("est_id=", signif(a$est, 8)),
                  paste0("est_full=", signif(b$est, 8)),
                  paste0("est_ref=", signif(ref$est, 8)),
                  paste0("sd_ic_hat=", signif(b$sd_idcond, 6))),
                collapse = "\t"),
          "\n", sep = "", file = out, append = TRUE)
    }
  }
  message("grid row ", g, " done")
}
