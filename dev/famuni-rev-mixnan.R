# Reviewer probe: the GEV mixture of dev/famuni-mixprobe.R (seed 35),
# fitted with the lane's soft support and with a NaN outside it.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(35); n <- 400
rg <- function(n, m, s, k) m + s * (rexp(n)^(-k) - 1) / k
d <- data.frame(yg = c(rg(n/2, 0, 1, 0.1), rg(n/2, 6, 1, 0.1)))
d$yf <- c(exp(0.5) / gamma(1 - 1/5) * rexp(n/2)^(-1/5),
          exp(2.5) / gamma(1 - 1/5) * rexp(n/2)^(-1/5))
d$yc <- c(rexp(n/2) + 1, rexp(n/2) * 3 + 8)
cat("range yg", range(d$yg), "\n")

nan_fam <- function() {
  fam <- gen_extreme_value()
  orig <- fam[["lpdf"]]
  fam[["lpdf"]] <- function(y, dpars, aterms) {
    g <- frmtmb:::gev_parts(y, dpars)
    t <- 1 + dpars[["xi"]] * (y - dpars[["mu"]]) / dpars[["sigma"]]
    # undo the lane's slope, and make an outside row NaN, as brms's is
    orig(y, dpars, aterms) + 1e8 * g[["out"]] + 0 * log(t)
  }
  fam
}
run <- function(lab, fam, v) {
  ws <- character(0)
  t0 <- proc.time()[[3]]
  fm <- tryCatch(withCallingHandlers(
    frm(as.formula(paste(v, "~ 1")), family = mixture(fam, fam), data = d),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fm, "error")) {
    cat(lab, v, "ERROR:", conditionMessage(fm), "\n"); return(invisible())
  }
  cat(sprintf("%s %s: logLik %.4f conv %d msg '%s' (%.1f s)\n", lab, v,
              as.numeric(logLik(fm)), fm$opt$convergence, fm$opt$message,
              proc.time()[[3]] - t0))
  cf <- fixef(fm, flatten = TRUE)
  print(signif(cf, 6))
  for (k in 1:2) {
    mu <- cf[[paste0("mu", k, "_(Intercept)")]]
    sg <- exp(cf[[paste0("sigma", k, "_(Intercept)")]])
    xi <- expm1(cf[[paste0("xi", k, "_(Intercept)")]])
    t <- 1 + xi * (d[[v]] - mu) / sg
    cat(sprintf("   component %d: end %s %.4f, rows outside %d, min t %.3e\n",
                k, if (xi > 0) "lower" else "upper", mu - sg / xi,
                sum(t <= 0), min(t)))
  }
  if (length(ws)) cat("   WARN:", unique(substr(ws, 1, 160)), sep = "\n   ")
  cat("\n")
  invisible(fm)
}
for (v in c("yg", "yc")) {
  run("soft", gen_extreme_value(), v)
  run("NaN ", nan_fam(), v)
}
