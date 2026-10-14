# Lane ci701: the boundary-correlation fit of test-review-v29.R:356
# (set.seed(61)), refitted with the data moved by a few ulps, to find
# the branch the macOS runner took.
#
#   Rscript dev/ci701-bcor-probe.R <lib|base> <K> data [k0]
#
# Fits k = k0..K: k = 0 is the test's data, k > 0 has every x, y1, y2
# times (1 + j 2^-52), j drawn per row from 0..4 by seed 1000 + k. The
# third argument is kept for the record and must be "data"; moved
# starts are dev/ci701-bcor-start.R.
args <- commandArgs(TRUE)
lib <- args[1]; K <- as.integer(args[2]); mode <- args[3]
k0 <- if (length(args) > 3) as.integer(args[4]) else 0L
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
cat("blas:", La_library(), extSoftVersion()[["BLAS"]], "\n")

set.seed(61)
ng <- 12
u <- stats::rnorm(ng, 0, 0.8)
per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
d0 <- dd

one <- function(k) {
  dd <- d0
  if (k > 0) {
    set.seed(1000 + k)
    for (v in c("x", "y1", "y2")) {
      dd[[v]] <- dd[[v]] * (1 + sample(0:4, nrow(dd), TRUE) * 2^-52)
    }
  }
  msgs <- character(0); warns <- character(0)
  fit <- withCallingHandlers(
    frm(mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
             bf(y2 ~ x + (1 | q | id)) + gaussian()), data = dd),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  sdr <- frmtmb:::sdr_of(fit)
  lost <- sdr$se_lost
  cw <- character(0)
  ci <- withCallingHandlers(confint_varcorr(fit), warning = function(w) {
    cw <<- c(cw, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  th <- fit$estimates$theta
  sdrow <- ci[ci$type == "sd", ]
  cr <- ci[ci$type == "cor", ]
  cat(sprintf(
    "RES k=%d code=%d msg=%s it=%d ll=%.10f th=%s lost=%s ok_sd=%s cor=%.8f cw=%d nmsg=%d nwarn=%d\n",
    k, fit$opt$convergence, gsub(" ", "_", fit$opt$message),
    as.integer(fit$opt$iterations %||% NA), as.numeric(logLik(fit)),
    paste(sprintf("%.6f", th), collapse = ","),
    if (length(lost)) paste(names(lost), lost, sep = ":", collapse = ",")
    else "none",
    all(is.finite(c(sdrow$lwr, sdrow$upr))), cr$estimate, length(cw),
    length(msgs), length(warns)))
  for (m in unique(c(msgs, warns, cw))) cat("   COND:", substr(m, 1, 300), "\n")
  invisible(fit)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
for (k in k0:K) one(k)
