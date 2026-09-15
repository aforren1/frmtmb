## Calibrating the suite test for item 2.6.
##
## The replicate sweep runs at the realistic design and takes half an
## hour a rung. The ordinary test suite cannot pay that, so the
## assertion that ships has to be a design small enough to run in
## seconds AND large enough that the width gap is not a coin toss.
##
## What this measures, over seeds, for a grid of designs:
##
##   ratio = width(cond + (1|id) + (1|id:cond)) / width(cond + (1|id))
##
## in a TREATMENT arm where the simulator has a subject-by-condition
## effect and in a NULL arm where it does not. The shipped assertion
## compares the two, so what matters is the gap between the arm minima
## and maxima, not either value.
##
## Run:
##   COH_SEEDS=20 Rscript dev/coh-calib.R > dev/coh-calib-log.txt
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
source("dev/coh-sim.R")

## No frequency axis: the contrast is what is being measured, so the
## rows inside a subject-by-condition cell are replicate draws and
## nothing in the truth depends on which one it is.
coh_flat <- function(seed, n_sub, n_rep, sd_idcond, sd_id = 0.35,
                     b0 = -0.6, b_cond = 0.5, nseg = 8L) {
  set.seed(seed)
  d <- expand.grid(rep = seq_len(n_rep), cond = factor(c("a", "b")),
                   id = factor(seq_len(n_sub)))
  u_id <- stats::rnorm(n_sub, 0, sd_id)
  u_ic <- stats::rnorm(n_sub * 2L, 0, sd_idcond)
  ic <- as.integer(d$id) + n_sub * (as.integer(d$cond) - 1L)
  eta <- b0 + b_cond * (d$cond == "b") +
    u_id[as.integer(d$id)] + u_ic[ic]
  w <- coupling_draw(rep(exp(0.3), nrow(d)), rep(exp(0.1), nrow(d)),
                     stats::plogis(eta), rep(0.4, nrow(d)), nseg)
  cbind(d, w)
}

seeds <- as.integer(Sys.getenv("COH_SEEDS", "20"))
grid <- expand.grid(n_sub = c(16L, 24L), n_rep = c(6L, 10L),
                    sd_ic = c(0.35, 0.5))
out <- "dev/coh-calib.tsv"
for (g in seq_len(nrow(grid))) {
  for (arm in c("treat", "null")) {
    sd_ic <- if (identical(arm, "null")) 0 else grid$sd_ic[g]
    for (s in seq_len(seeds)) {
      d <- coh_flat(4200000L + 1000L * g + s, grid$n_sub[g],
                    grid$n_rep[g], sd_ic)
      t0 <- Sys.time()
      a <- try(coh_fit_one("cond + (1 | id)", d), silent = TRUE)
      b <- try(coh_fit_one("cond + (1 | id) + (1 | id:cond)", d),
               silent = TRUE)
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (inherits(a, "try-error") || inherits(b, "try-error")) {
        cat(sprintf("g=%d arm=%s seed=%d ERROR\n", g, arm, s),
            file = out, append = TRUE)
        next
      }
      cat(paste(c(paste0("g=", g), paste0("n_sub=", grid$n_sub[g]),
                  paste0("n_rep=", grid$n_rep[g]),
                  paste0("sd_ic_true=", sd_ic),
                  paste0("arm=", arm), paste0("s=", s),
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
  }
  message("grid row ", g, " done")
}
