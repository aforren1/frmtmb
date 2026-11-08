# Lane latperf: time, memory and saved size of the readers of a latent
# residual's fitted rows, one design per process.
#
# Usage: Rscript dev/latperf-prof.R <arm> <n> <design> [dense] [rprof]
#   arm     base (rellib-r18 alone) or lane (wt-latperf-lib first)
#   n       rows
#   design  one: y ~ x + ar(t, g, cov = TRUE), one series, poisson
#           two: the same plus (1 | id) with 50 levels, and fitted() on 50
#                new rows
#   dense   "dense" moves frmtmb.latres_se_rows and frmtmb.latres_dense_rows
#           to Inf, so the base build takes its dense route at any n
#   rprof   "rprof" writes an Rprof summary of fitted()
# Seed 1, the generator of dev/review-latres-fittedprof.R.
args <- commandArgs(TRUE)
arm <- args[1]
n <- as.integer(args[2])
design <- args[3]
dense <- "dense" %in% args
rprof <- "rprof" %in% args
base <- c("C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lane") c("C:/Users/adf44/source/r/wt-latperf-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
if (dense) {
  options(frmtmb.latres_se_rows = Inf, frmtmb.latres_dense_rows = Inf)
}
peak_ws <- function() {
  out <- system2("powershell", c("-NoProfile", "-Command",
    sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid())),
    stdout = TRUE)
  as.numeric(out[length(out)]) / 2^20
}
heap_mb <- function() {
  g <- gc()
  # "max used" in Mb of Ncells and Vcells since the last reset
  sum(g[, ncol(g)])
}
cat(sprintf("lib: %s frmtmb %s\n", dirname(dirname(system.file(package =
  "frmtmb"))), format(packageVersion("frmtmb"))))
cat(sprintf("arm %s n %d design %s dense %s\n", arm, n, design, dense))
set.seed(1)
d <- data.frame(t = seq_len(n), g = factor(1), x = rnorm(n),
                id = factor(sample.int(50, n, TRUE)))
e <- as.numeric(arima.sim(list(ar = .5), n, sd = .5))
u <- rnorm(50, 0, 0.3)
d$y <- rpois(n, exp(.5 + .3 * d$x + e +
                     if (design == "two") u[d$id] else 0))
f <- if (design == "two") {
  y ~ x + (1 | id) + ar(t, g, cov = TRUE)
} else {
  y ~ x + ar(t, g, cov = TRUE)
}
invisible(gc(reset = TRUE))
tf <- system.time(fit <- suppressWarnings(
  frm(f, family = poisson(), data = d)))[["elapsed"]]
cat(sprintf("RESULT step fit time %.2f heap %.0f ws %.0f\n", tf, heap_mb(),
            peak_ws()))
cat(sprintf("RESULT step saved_size_before MB %.2f\n",
            length(serialize(fit, NULL)) / 2^20))
step <- function(name, expr) {
  invisible(gc(reset = TRUE))
  msg <- character(0)
  t <- system.time(v <- tryCatch(withCallingHandlers(expr,
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      msg <<- c(msg, paste("ERROR:", conditionMessage(e)))
      NULL
    }))[["elapsed"]]
  se <- if (is.list(v)) v$se.fit else if (is.matrix(v)) v[, "Est.Error"]
  cat(sprintf(paste("RESULT step %s time %.2f heap %.0f ws %.0f se_na %d",
                    "se_first %.10g se_med %.10g\n"),
              name, t, heap_mb(), peak_ws(), sum(is.na(se)), se[1],
              stats::median(se, na.rm = TRUE)))
  if (length(msg)) cat("  message:", substr(msg[1], 1, 160), "\n")
  invisible(v)
}
if (rprof) {
  tmp <- tempfile()
  Rprof(tmp, interval = 0.02)
}
a <- step("fitted", suppressWarnings(fitted(fit)))
if (rprof) {
  Rprof(NULL)
  s <- summaryRprof(tmp)
  print(head(s$by.self, 10))
  bt <- s$by.total
  print(head(bt[order(-bt$total.time), c("total.time", "self.time")], 30))
}
step("residuals", suppressWarnings(residuals(fit, summary = TRUE)))
step("linpred_se", frm_linpred(fit, se.fit = TRUE))
raw <- serialize(fit, NULL)
cat(sprintf("RESULT step saved_size MB %.2f cache %s\n", length(raw) / 2^20,
            paste(ls(fit$cache), collapse = ",")))
if (design == "two") {
  set.seed(2)
  nd <- data.frame(t = n + seq_len(50), g = factor(1), x = rnorm(50),
                   id = factor(sample.int(50, 50, TRUE), levels = 1:50))
  step("fitted_newdata", suppressWarnings(fitted(fit, newdata = nd)))
  step("linpred_newdata_se", frm_linpred(fit, newdata = nd, se.fit = TRUE))
}
cat(sprintf("RESULT step end ws %.0f\n", peak_ws()))
