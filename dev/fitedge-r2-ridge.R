# Fix round, m2: the ridge pass keyed on a stop warning (shipped) against the
# pass on any held warning (the reviewed build: nl_stop_warning() patched to
# `length(held) > 0`), in one process, order alternating by row, on the
# review's `mvs` and `mv4` groups. The data and models are
# dev/review-fitedge-nl.R's own functions, read from that file:
# gen(seed * 100 + 7, slope), seeds 1:4, slopes +-0.6. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ex <- parse(file.path(FE_WT, "dev/review-fitedge-nl.R"))
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("gen", "pars_of", "add_lf", "one_body",
                                    "one_fam", "models", "fitm")) {
    eval(e, globalenv())
  }
}
ns <- asNamespace("frmtmb")
orig <- get("nl_stop_warning", envir = ns)
use <- function(new) {
  unlockBinding("nl_stop_warning", ns)
  assign("nl_stop_warning", if (new) orig else function(held) length(held) > 0L,
         envir = ns)
  lockBinding("nl_stop_warning", ns)
}
one <- function(new, m, d) {
  use(new)
  t0 <- proc.time()[[3]]
  r <- fitm(m, d)
  list(t = proc.time()[[3]] - t0,
       ll = if (inherits(r$fit, "error")) NA else as.numeric(logLik(r$fit)),
       w = sort(substr(r$warn, 1, 60)))
}
i <- 0L
for (grp in c("mvs", "mv4")) {
  m <- models(grp)[[1]]
  tn <- to <- 0
  same <- 0L
  wsame <- 0L
  for (seed in 1:4) for (slope in c(0.6, -0.6)) {
    i <- i + 1L
    d <- gen(seed * 100 + 7, slope)
    if (i %% 2L) {
      a <- one(TRUE, m, d); b <- one(FALSE, m, d)
    } else {
      b <- one(FALSE, m, d); a <- one(TRUE, m, d)
    }
    tn <- tn + a$t
    to <- to + b$t
    same <- same + isTRUE(abs(a$ll - b$ll) <= 1e-8 * abs(b$ll))
    wsame <- wsame + identical(a$w, b$w)
    cat(sprintf("RG %s seed %d slope %.1f stop-keyed %.6f (%.1f s) any-warning %.6f (%.1f s)\n",
                grp, seed, slope, a$ll, a$t, b$ll, b$t))
  }
  cat(sprintf("RGSUM %s rows 8: stop-keyed %.1f s, any-warning %.1f s, ratio %.2f; logLik equal %d of 8, warnings equal %d of 8\n",
              grp, tn, to, tn / to, same, wsame))
}
cat("DONE\n")
