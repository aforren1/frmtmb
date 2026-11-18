# Lane squant: does quantile_elf() with an observation-level intercept
# and weights above 1 run away, as the families of
# unbounded_spike_dpar (R/weights-unbounded.R) do? Weights 0.9 and 1.2,
# the construction of dev/fuzzfix-spike-families.R; seeds 1 to 3.
# The Laplace objective along log sigma, the rest re-optimized.
#
# Run: Rscript dev/squant-weights.R > dev/squant-out/weights.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
for (s in 1:3) {
  set.seed(s)
  n <- 60
  d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  for (w in c(0.9, 1.2)) {
    d$w <- w
    ws <- character()
    fit <- tryCatch(withCallingHandlers(
      frm(bf(y | weights(w) ~ x + (1 | obs), quantile = 0.5),
          family = quantile_elf(), data = d),
      warning = function(e) {
        ws <<- c(ws, substr(conditionMessage(e), 1, 90))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (inherits(fit, "error")) {
      cat(sprintf("seed %d w %.1f: ERROR %s\n", s, w,
                  substr(gsub("\\s+", " ", conditionMessage(fit)), 1, 160)))
      next
    }
    par <- fit$opt$par
    j <- which(names(par) == "betad")[1L]
    prof <- vapply(c(0, -2, -5, -10, -20), function(k) {
      o <- nlminb(par[-j], function(p) {
        full <- par
        full[-j] <- p
        full[j] <- par[j] + k
        v <- fit$obj$fn(full)
        if (is.finite(v)) v else 1e10
      })
      -o$objective
    }, 0)
    cat(sprintf(paste0("seed %d w %.1f: conv %d logLik %.3f log sigma %.3f |",
                       " profile at log sigma +0 -2 -5 -10 -20: %s | %s\n"),
                s, w, fit$opt$convergence, as.numeric(logLik(fit)), par[j],
                paste(sprintf("%.2f", prof), collapse = " "),
                paste(unique(ws), collapse = " | ")))
  }
}
