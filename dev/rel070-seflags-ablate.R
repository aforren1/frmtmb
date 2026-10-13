# rel070 seflags: attribution of the silent shortfalls of
# rel070-seflags-shortscan.R. In one process, each fit twice: as
# installed, and with optrap's response units switched off for this
# process only (predictor_units() rebound in the loaded namespace to
# return no unit). Nothing on disk changes.
#   Rscript rel070-seflags-ablate.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "| R", R.home(), "\n")
pu <- get("predictor_units", ns)
set_units <- function(on) {
  unlockBinding("predictor_units", ns)
  assign("predictor_units",
         if (on) pu else function(frame) numeric(0), envir = ns)
  lockBinding("predictor_units", ns)
}
cases <- list(c(6, 13), c(6, 15), c(7, 3), c(7, 9), c(7, 14), c(7, 17),
              c(7, 18), c(7, 20), c(8, 3))
for (cs in cases) {
  deg <- cs[1]; s <- cs[2]
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs, "+ (1 | g)"))
  m <- suppressWarnings(suppressMessages(lme4::lmer(
    stats::as.formula(paste("y ~ poly(x,", deg, ") + (1 | g)")), data = d,
    REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 1e5)))))
  out <- character()
  for (on in c(TRUE, FALSE)) {
    set_units(on)
    w <- character()
    f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
    l <- suppressWarnings(ns$sdr_of(f)$se_lost)
    out <- c(out, sprintf("units %-3s: code %d short %.4f warnings %d lost %d",
                          if (on) "on" else "off", f$opt$convergence,
                          as.numeric(logLik(m)) - as.numeric(logLik(f)),
                          length(w), length(l)))
  }
  set_units(TRUE)
  cat(sprintf("degree %d seed %2d | unit %s | %s | %s\n", deg, s,
              format(pu(f$frame)), out[1], out[2]))
}
