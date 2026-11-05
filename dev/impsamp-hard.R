# The designs the old guard and the old tests called unidentified or
# stalled, and very small draw counts, on the fixed build. Usage:
#   Rscript dev/impsamp-hard.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("dev/impsamp-agh.R")
show <- function(lab, form, fam, d, Ns, shape = NULL, famname = NULL,
                 ctl = list()) {
  cat("\n---", lab, "---\n")
  lap <- suppressWarnings(suppressMessages(frm(form, family = fam, data = d)))
  cat(sprintf("Laplace logLik %.5f convergence %d par %s\n",
              as.numeric(logLik(lap)), lap$opt$convergence,
              paste(format(lap$opt$par, digits = 4), collapse = " ")))
  ref <- NULL
  if (!is.null(shape)) {
    ref <- impsamp_ref(famname, shape, d, lap)
    cat(sprintf("reference (own AGH, optim) logLik %.5f par %s\n", ref$ll,
                paste(format(ref$par, digits = 4), collapse = " ")))
  }
  for (N in Ns) {
    ws <- character(0)
    r <- tryCatch(withCallingHandlers(suppressMessages(frm(
      form, family = fam, data = d, importance = N,
      control = do.call(frmtmb_control, ctl))),
      warning = function(w) {
        ws <<- c(ws, substr(conditionMessage(w), 1, 110))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (inherits(r, "error")) {
      cat(sprintf("N=%d ERROR %s\n", N, substr(conditionMessage(r), 1, 400)))
    } else {
      im <- r$importance
      ex <- if (!is.null(shape)) {
        tryCatch(impsamp_exact_at(famname, shape, d, r),
                 error = function(e) NA_real_)
      } else NA_real_
      cat(sprintf(paste0("N=%d ll %.5f exact-at-own %.5f mcse %.4f essmin %.3f",
                         " rounds %d capped %s moves %s par %s\n"),
                  N, as.numeric(logLik(r)), ex, im$mcse, im$ess_min, im$rounds,
                  im$capped, paste(format(im$moves, digits = 3), collapse = ","),
                  paste(format(r$opt$par, digits = 4), collapse = " ")))
      if (length(ws)) cat("   warn:", ws, sep = " | ", "\n")
    }
  }
}
Ns <- c(4L, 8L, 16L, 32L, 64L, 200L, 500L, 2000L)
# test-importance.R's "diverging" design, seed 11
set.seed(11)
ng <- 40L
per <- 4L
g <- rep(seq_len(ng), each = per)
x <- rnorm(ng * per)
b <- cbind(rnorm(ng, 0, 1), rnorm(ng, 0, 0.7))
dd <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.8 * x + b[g, 1] +
                                                  b[g, 2] * x)),
                 x = x, g = factor(g))
dd$ntr <- 1L
show("test-importance diverging design", bf(y ~ x + (x | g)), bernoulli(), dd,
     Ns, "slope", "bernoulli")
# test-importance.R's stall design: sd = 0, 12 groups of 3
set.seed(1)
ng <- 12L
per <- 3L
g <- factor(rep(seq_len(ng), each = per))
x <- rnorm(ng * per)
ds <- data.frame(y = rbinom(ng * per, 1, plogis(-0.3 + 0.7 * x)), x = x, g = g,
                 ntr = 1L)
show("stall design", bf(y ~ x + (1 | g)), bernoulli(), ds, Ns, "int",
     "bernoulli")
# the doc's hard design: forty groups of three binary rows, sd 2
for (sd_seed in c(1L, 3L, 7L)) {
  set.seed(sd_seed)
  ng <- 40L
  per <- 3L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 2)
  dh <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.5 * x + u[g])),
                   x = x, g = g, ntr = 1L)
  show(paste("hard design, sd 2, seed", sd_seed), bf(y ~ x + (1 | g)),
       bernoulli(), dh, Ns, "int", "bernoulli")
}
# an unidentified covariance: q = 3 on 2 rows per group
set.seed(5)
ng <- 30L
g <- factor(rep(seq_len(ng), each = 2L))
x <- rnorm(2L * ng)
du <- data.frame(y = rbinom(2L * ng, 1, 0.4), x = x, g = g, ntr = 1L)
show("unidentified, null data, 30 groups of 2", bf(y ~ x + (x | g)),
     bernoulli(), du, Ns, "slope", "bernoulli")
