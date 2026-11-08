# Probe: where marginaleffects calls get_coef() and get_vcov(), by caller,
# and how long one get_predict() takes on a mixed fit (data seed 7411).
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
set.seed(7411)
J <- 30; m <- 40; n <- J * m
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(J), each = m)))
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g]))
fit <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = d))
ns <- asNamespace("frmtmb")
log <- new.env()
log$calls <- character(0)
wrap <- function(nm) {
  f0 <- get(paste0(nm, ".frmtmb_fit"), envir = ns)
  function(model, ...) {
    sc <- sys.calls()
    callers <- vapply(sc, function(cl) {
      f <- cl[[1L]]
      if (is.name(f)) as.character(f) else paste(deparse(f), collapse = "")
    }, "")
    callers <- callers[!callers %in% c("eval", "do.call", "do_call",
                                       "withCallingHandlers", "tryCatch",
                                       "tryCatchList", "tryCatchOne",
                                       "doTryCatch", "suppressWarnings",
                                       "suppressMessages", "quiet", "run",
                                       "FUN", "lapply", "vapply", "sapply")]
    log$calls <- c(log$calls, paste(nm, "[", paste(names(list(...)),
                                                   collapse = ","), "] <-",
                                    paste(tail(callers, 4), collapse = " < ")))
    f0(model, ...)
  }
}
for (nm in c("get_coef", "get_vcov")) {
  registerS3method(nm, "frmtmb_fit", wrap(nm),
                   envir = asNamespace("marginaleffects"))
}
run <- function(lab, e) {
  log$calls <- character(0)
  r <- quiet(e)
  cat("\n==", lab, "\n")
  print(table(log$calls))
}
run("predictions", predictions(fit, newdata = d[1:3, ]))
run("avg_predictions re_formula NA", avg_predictions(fit, re_formula = NA))
run("avg_slopes", avg_slopes(fit, variables = "x"))
run("avg_comparisons", avg_comparisons(fit, variables = "x"))
run("hypotheses", hypotheses(fit))
run("hypotheses b2", hypotheses(fit, hypothesis = "b2 = 0"))
run("predictions hyp", predictions(fit, newdata = d[1:3, ],
                                   hypothesis = "b1 = b2"))
cat("\ntiming one get_predict on", n, "rows:\n")
t0 <- proc.time()[["elapsed"]]
k <- 0L
while (proc.time()[["elapsed"]] - t0 < 3) {
  get_predict(fit, newdata = d)
  k <- k + 1L
}
cat(sprintf("%.4f s per call over %d calls\n",
            (proc.time()[["elapsed"]] - t0) / k, k))
t0 <- proc.time()[["elapsed"]]
invisible(quiet(avg_slopes(fit, variables = "x")))
cat(sprintf("avg_slopes: %.2f s\n", proc.time()[["elapsed"]] - t0))
t0 <- proc.time()[["elapsed"]]
invisible(quiet(avg_predictions(fit)))
cat(sprintf("avg_predictions: %.2f s\n", proc.time()[["elapsed"]] - t0))
