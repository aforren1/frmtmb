# Lane famend, 3060: frechet() fits near nu = 1. For each fit, nu - 1,
# the optimizer code, whether the SEs of the mean intercept and nu are
# finite, and d logLik / d nu at the stop (the GEV check's slope, read
# off the objective's gradient carried to the nu scale). Edge fits and
# interior maxima, to place a threshold between them.
#   Rscript dev/famend-probe-frechet.R <lib or "base"> <design> [seeds]
# designs: n300 (nu0 = 1.1, n = 300, dev/famuni-rev-frechet11.R),
# near (nu0 = 1.02, n = 3000), far (nu0 = 1.5, n = 300), heavy (nu0 =
# 1.05, n = 20000), mult (nu0 = 1.1, n = 300, nu ~ x2)
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
des <- a[2]
seeds <- if (length(a) > 2) eval(parse(text = a[3])) else 1:60
cat("lib", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "design", des, "\n")
cfg <- switch(des, n300 = c(1.1, 300), near = c(1.02, 3000),
              far = c(1.5, 300), heavy = c(1.05, 20000),
              mult = c(1.1, 300))
for (s in seeds) {
  set.seed(9000 + s)
  nu0 <- cfg[1]
  n <- cfg[2]
  d <- data.frame(x = rnorm(n), x2 = rnorm(n))
  d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(n)^(-1 / nu0)
  ws <- character(0)
  ms <- character(0)
  fo <- if (des == "mult") bf(y ~ x, nu ~ x2) else bf(y ~ x)
  f <- tryCatch(withCallingHandlers(frm(fo, family = frechet(), data = d),
    warning = function(e) {
      ws <<- c(ws, conditionMessage(e)); invokeRestart("muffleWarning")
    },
    message = function(e) {
      ms <<- c(ms, conditionMessage(e)); invokeRestart("muffleMessage")
    }), error = function(e) conditionMessage(e))
  if (is.character(f)) {
    cat(sprintf("seed %d: ERROR %s\n", 9000 + s, substr(f, 1, 100)))
    next
  }
  nu <- frm_linpred(f, dpar = "nu", type = "response")
  lp <- f$frame[["linpreds"]][[ns$linpred_key(f$spec$responses[[1]]$resp_name,
                                               "nu")]]
  dll <- if (des != "mult") ns$gev_xi_slope(f, lp) else NA
  gain <- if (des != "mult" && exists("frechet_ridge_gain", ns)) {
    ns$frechet_ridge_gain(f, f$spec$responses[[1]]$resp_name, lp)
  } else NA
  dll <- paste0(format(dll, digits = 3), " gain ", format(gain, digits = 3))
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  g <- f$obj$gr(f$opt$par)
  cat(sprintf(paste0("seed %d: code %d min(nu - 1) %.3g dll %s maxgrad %.2g ",
                     "se finite %d of %d | W:%s | M:%s\n"),
              9000 + s, f$opt$convergence, min(nu) - 1,
              format(dll, digits = 3), max(abs(g)),
              sum(is.finite(se)), length(se),
              if (length(ws)) substr(gsub("\\s+", " ", ws[1]), 1, 70) else "",
              if (length(ms)) substr(gsub("\\s+", " ", ms[1]), 1, 70) else ""))
}
cat("done\n")
