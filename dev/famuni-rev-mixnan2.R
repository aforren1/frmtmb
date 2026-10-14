# Reviewer probe: the GEV mixture of dev/famuni-mixprobe.R (seed 35),
# fitted with the lane's soft support and, by replacing gev_parts() in
# the namespace, with a NaN outside the support (brms's behavior).
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

soft_parts <- frmtmb:::gev_parts
nan_parts <- function(q, dpars) {
  x <- (q - dpars[["mu"]]) / dpars[["sigma"]]
  z <- dpars[["xi"]] * x
  out <- log1p(z) - log1p(z)           # 0 inside, NaN outside; 0 * x folds to 0 on the tape
  zc <- ad_max(z + out, -1 + 1e-15)
  list(l = log1p(zc), h = x * log1p_ratio(zc), out = out)
}
environment(nan_parts) <- asNamespace("frmtmb")

run <- function(lab, v, fam = gen_extreme_value()) {
  ws <- character(0)
  t0 <- proc.time()[[3]]
  fm <- tryCatch(withCallingHandlers(
    frm(as.formula(paste(v, "~ 1")), family = mixture(fam, fam), data = d),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fm, "error")) {
    cat(lab, v, "ERROR:", conditionMessage(fm), "\n\n"); return(invisible())
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
  assignInNamespace("gev_parts", soft_parts, "frmtmb")
  run("soft", v)
  assignInNamespace("gev_parts", nan_parts, "frmtmb")
  # the replacement must reach the density: a row outside is NaN
  chk <- gen_extreme_value()$lpdf(c(0, -100), list(mu = 0, sigma = 1,
                                                   xi = 0.5), list())
  cat("NaN check, lpdf at an inside and an outside row:", chk, "\n")
  run("NaN ", v)
  run1 <- tryCatch(frm(as.formula(paste(v, "~ 1")),
                       family = gen_extreme_value(), data = d),
                   error = function(e) e)
  if (!inherits(run1, "error")) {
    cat("single GEV under NaN:", as.numeric(logLik(run1)),
        run1$opt$convergence, "\n")
  } else cat("single GEV under NaN: ERROR", conditionMessage(run1), "\n")
  assignInNamespace("gev_parts", soft_parts, "frmtmb")
  run1 <- frm(as.formula(paste(v, "~ 1")), family = gen_extreme_value(),
              data = d)
  cat("single GEV soft:     ", as.numeric(logLik(run1)),
      run1$opt$convergence, "\n\n")
}
