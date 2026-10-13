# Consolidation measurement, not a change to the package: the release
# build's newton_polish() accepts a step whenever the objective
# decreases, so on negative binomial data at the Poisson limit (the
# shape's MLE infinite) a step predicted to gain more than 1e-6 that
# gains 3e-8 runs the log shape from 20 to 67, into an indefinite
# Hessian, and every SE is lost behind a large-gradient warning
# (dev/test-backlog.md, "Found at the consolidation"). This patches, in
# this process only, a candidate gate: a step is kept only when its
# actual gain is at least `GATE` (default 0.25) of the predicted one.
#
#   Rscript dev/rel070-polishgate.R negbin <arm> <seeds expr>
#   Rscript dev/rel070-polishgate.R poly <arm> <label>   # polystudy
# <arm>: "asis" (the release build) or "gate". POLISHGATE_LIB=<lib>
# runs another build (rellib-r7, 0.69.1, for the before arm).
a <- commandArgs(TRUE)
mode <- a[1]
arm <- a[2]
.libPaths(c(Sys.getenv("POLISHGATE_LIB", "C:/Users/adf44/source/r/rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
gate <- as.numeric(Sys.getenv("GATE", "0.25"))
if (identical(arm, "gate")) {
  # the function's source in the release tree, which the installed
  # build was made from, with the one line added
  src <- readLines("R/fit.R", warn = FALSE)
  s0 <- grep("^newton_polish <- function", src)
  e0 <- s0 - 1L + which(src[s0:length(src)] == "}")[1]
  b <- src[s0:e0]
  i <- grep("if (!is.finite(vq) || vq >= v) break", b, fixed = TRUE)
  stopifnot(length(s0) == 1L, length(i) == 1L)
  b <- append(b, sprintf("    if (v - vq < %s * gain) break", gate),
              after = i)
  f <- eval(parse(text = b)[[1]][[3]], envir = ns)
  utils::assignInNamespace("newton_polish", f, "frmtmb")
  stopifnot(any(grepl("v - vq <", deparse(body(
    get("newton_polish", asNamespace("frmtmb")))), fixed = TRUE)))
}
cat("mode", mode, "arm", arm, "gate", gate, "lib",
    dirname(find.package("frmtmb")), "\n")
if (mode == "negbin") {
  seeds <- eval(parse(text = a[3]))
  out <- NULL
  for (s in seeds) {
    # dev/seflags-rev3-falarm.R's negbin_large_shape data, its RNG order
    set.seed(s)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    set.seed(s)
    n <- 150
    x <- rnorm(n)
    yy <- 0.2 + 0.5 * x + rnorm(n)
    mu <- plogis(-2 + 0.5 * x)
    yb <- rbeta(n, mu * 200, (1 - mu) * 200)
    yb[sample(n, 3)] <- 0
    yb2 <- rbeta(n, 0.98 * 20, 0.02 * 20)
    yb2[sample(n, 2)] <- 1
    yt <- 1 + 0.5 * x + rt(n, df = 30)
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
      seed = s, code = fit$opt$convergence,
      ll = as.numeric(logLik(fit)), shape = fit$opt$par[3],
      se_int = se[1], se_x = se[2], nwarn = length(w),
      steps = if (is.null(pol)) 0 else pol[["steps"]],
      gain = if (is.null(pol)) 0 else pol[["gain"]]))
  }
  print(out, row.names = FALSE, digits = 7)
  cat("SUMMARY", arm, ": fits", nrow(out), "| warned", sum(out$nwarn > 0),
      "| NaN se of b", sum(!is.finite(out$se_int) | !is.finite(out$se_x)),
      "| polished", sum(out$steps > 0), "\n")
} else {
  # dev/optrap-p1-polystudy.R reads commandArgs(TRUE)[1:2] as <lib>
  # <label>; it sees this script's arguments, so its library line names
  # the arm and the label
  source("dev/optrap-p1-polystudy.R", echo = FALSE)
}
