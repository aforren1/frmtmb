# Lane seflags, punch round 1, B1 under the reference BLAS: on
# `y ~ x + s(x) + (1 | g)` (dev/seflags-alias.R's data, mean 1e-5, and
# the review's alias plus an identified pair z2 = z1 + N(0, 1e-4)), the
# hypothesis b_x + t * bs_sx_1 is estimable only at t = b, the slope of
# the smooth's fixed column on x. Sweep t = b (1 - delta); print the
# Est.Error and the warnings. With "prepunch", the lane's se-check.R as
# it was before the punch round is evaluated over the namespace.
#   Rscript dev/seflags-b1ref.R [prepunch]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
if (length(a) && a[1] == "prepunch") {
  tmp <- new.env()
  sys.source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-out/prepunch-se-check.R", envir = tmp)
  for (nm in ls(tmp, all.names = TRUE)) if (exists(nm, envir = ns, inherits = FALSE)) {
    unlockBinding(nm, ns); f <- get(nm, envir = tmp)
    if (is.function(f)) environment(f) <- ns
    assign(nm, f, envir = ns)
  }
}
cat("lib", find.package("frmtmb"), "| arm", if (length(a)) a[1] else "lane", "| R", R.home(), "\n")
mk <- function(pair) {
  set.seed(1)
  n <- 200
  x <- stats::rnorm(n)
  x <- x - mean(x) + 1e-5
  d <- data.frame(x = x, g = factor(rep(1:10, n / 10)))
  d$z1 <- stats::rnorm(n)
  d$z2 <- d$z1 + stats::rnorm(n, 0, 1e-4)
  d$y <- 1 + 0.5 * d$x + sin(2 * d$x) + stats::rnorm(10, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.3) + if (pair) 0.3 * (d$z1 + d$z2) else 0
  d
}
for (pair in c(FALSE, TRUE)) {
  d <- mk(pair)
  fo <- if (pair) y ~ x + s(x) + z1 + z2 + (1 | g) else y ~ x + s(x) + (1 | g)
  fit <- suppressWarnings(frm(fo, data = d))
  X <- fit$frame$linpreds[[1]]$X
  b <- unname(coef(stats::lm(as.numeric(X[, "s(x).fx1"]) ~ d$x))[2])
  sdr <- suppressWarnings(ns$sdr_of(fit))
  cat(sprintf("\npair %s | b %.10g | lost %s | flat tol %s\n", pair, b,
              paste(names(sdr$se_lost), collapse = ","),
              format(sdr$se_null_flat$tol)))
  for (dl in c(0, 1e-4, 1e-3, 1e-2, 3e-2, 0.1)) {
    t <- b * (1 - dl)
    w <- character()
    h <- withCallingHandlers(hypothesis(fit, sprintf("b_x + %.15g * bs_sx_1 = 0", t), class = NULL),
      warning = function(z) { w <<- c(w, conditionMessage(z)); invokeRestart("muffleWarning") })
    cat(sprintf("  delta %-6g | Est.Error %-10s | warnings %d\n", dl,
                format(h$hypothesis$Est.Error, digits = 6), length(w)))
  }
}
