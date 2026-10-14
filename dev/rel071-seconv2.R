# Candidate starts for test-se-conv.R's sd-near-zero case on one build:
# which reproduce a coded stop short of the maximum whatever the
# default start of mu and sigma is? `start = list(theta = -15)` alone
# lets lane optstart's least-squares start place beta and sigma, so the
# fit converges at code 0 to the maximum given sd = 0.
#   Rscript dev/rel071-seconv2.R <lib> [label]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("==", if (length(a) > 1) a[2] else "", "frmtmb",
    format(packageVersion("frmtmb")), dirname(find.package("frmtmb")), "\n")
conv_phrase <- "Optimizer did not report convergence"
run <- function(expr) {
  ws <- character(0)
  v <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = ws)
}
starts <- list(
  "theta" = list(theta = -15),
  "beta 0, theta" = list(beta = c(0, 0), theta = -15),
  "beta 0, betad 0, theta" = list(beta = c(0, 0), betad = 0, theta = -15))
for (nm in names(starts)) for (s in 1:4) {
  set.seed(s)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.4)
  ref <- run(frm(y ~ x + (1 | g), data = d))$v
  r <- run(frm(y ~ x + (1 | g), data = d, start = starts[[nm]],
               control = frmtmb_control(optCtrl = list(iter.max = 3,
                                                       eval.max = 1000))))
  f <- r$v
  w <- grep(conv_phrase, r$w, fixed = TRUE, value = TRUE)
  cat(sprintf(paste0("%-24s seed %d: code %d, gap %9.4f, conv warning %d,",
                     " 'short of the maximum' %s, 'at the maximum' %s,",
                     " 'may not have converged' %d\n"),
              nm, s, f$opt$convergence,
              as.numeric(logLik(f)) - as.numeric(logLik(ref)), length(w),
              any(grepl("short of the maximum", w, fixed = TRUE)),
              any(grepl("so the estimates are at the maximum", w, fixed = TRUE)),
              length(grep("The fit may not have converged", r$w, fixed = TRUE))))
}
