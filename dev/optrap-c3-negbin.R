# Lane optrap, consolidation fixes: the 0.70.0 consolidation's negative
# binomial check (its dev/rel070-polishgate.R negbin mode, the data of
# dev/seflags-rev3-falarm.R's negbin_large_shape in their RNG order),
# 200 seeds: fits that warn, that lose a coefficient SE, that the Newton
# steps move, and the largest gain.
#   Rscript dev/optrap-c3-negbin.R <lib> <label> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
seeds <- if (length(a) > 2) eval(parse(text = a[3])) else 1:200
cat("lib", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "label", a[2], "\n")
out <- NULL
for (s in seeds) {
  set.seed(s)
  n <- 150
  x <- rnorm(n)
  invisible(rnorm(n))
  mu <- plogis(-2 + 0.5 * x)
  invisible(rbeta(n, mu * 200, (1 - mu) * 200)); invisible(sample(n, 3))
  invisible(rbeta(n, 0.98 * 20, 0.02 * 20)); invisible(sample(n, 2))
  invisible(rt(n, df = 30))
  dn <- data.frame(x = x, c = rpois(n, exp(1 + 0.3 * x)))
  w <- character()
  fit <- withCallingHandlers(
    frm(c ~ x, family = negbinomial(), data = dn),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage"))
  se <- suppressWarnings(fixef(fit)[, "Est.Error"])
  pol <- fit$opt$polish
  out <- rbind(out, data.frame(
    seed = s, code = fit$opt$convergence, ll = as.numeric(logLik(fit)),
    shape = fit$opt$par[3], se_int = se[1], se_x = se[2],
    nwarn = length(w), steps = if (is.null(pol)) 0 else pol[["steps"]],
    gain = if (is.null(pol)) 0 else pol[["gain"]]))
}
cat("SUMMARY", a[2], ": fits", nrow(out), "| warned", sum(out$nwarn > 0),
    "| NaN se of b", sum(!is.finite(out$se_int) | !is.finite(out$se_x)),
    "| polished", sum(out$steps > 0), "| largest gain",
    format(max(out$gain), digits = 3), "\n")
print(out[out$steps > 0 | !is.finite(out$se_int), ], row.names = FALSE,
      digits = 6)
