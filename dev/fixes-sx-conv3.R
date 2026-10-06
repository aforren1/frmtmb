# Lane fixes, punch round: how often a smooth fit ends without
# convergence or with non-finite standard errors, base vs lane, over
# gamSim seeds 1 to 20 and seven smooth formulas (140 fits per build).
#   Rscript dev/fixes-sx-conv3.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
forms <- list(y ~ s(x1) + s(x2), y ~ s(x1, bs = "cr", k = 6),
              y ~ s(x1, by = g) + g, y ~ s(x1, by = z), y ~ s(x1, x2),
              y ~ t2(x1, x2), y ~ s(x0) + s(x1) + s(x2) + s(x3))
rows <- list()
for (seed in 1:20) {
  set.seed(seed)
  d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, scale = 2,
                                     verbose = FALSE))
  d$z <- runif(200)
  d$g <- factor(sample(c("a", "b", "c"), 200, TRUE))
  for (k in seq_along(forms)) {
    nw <- 0L
    fit <- withCallingHandlers(frm(bf(forms[[k]]), data = d),
                               warning = function(x) {
                                 nw <<- nw + 1L
                                 invokeRestart("muffleWarning")
                               })
    se <- suppressWarnings(fixef(fit)[, "Est.Error"])
    rows[[length(rows) + 1L]] <- data.frame(
      seed = seed, form = k, conv = fit$opt$convergence,
      se_ok = all(is.finite(se)), warns = nw,
      ll = as.numeric(logLik(fit)))
  }
}
r <- do.call(rbind, rows)
saveRDS(r, sub("[.]R$", paste0("-", basename(dirname(find.package("frmtmb"))),
                             ".rds"), "dev/fixes-log/sx-conv3.R"))
cat(sprintf("fits %d; conv != 0: %d; non-finite fixef SE: %d; with a warning: %d\n",
            nrow(r), sum(r$conv != 0), sum(!r$se_ok), sum(r$warns > 0)))
print(r[r$conv != 0 | !r$se_ok, ])
