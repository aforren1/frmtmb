# Item 3477: a quadrature fit with a variance at zero. The one-block case
# of dev/fuzzfix-quad-boundary.R (seed 3, bernoulli, 20 groups of 6),
# under Laplace and under quadrature: messages, the sd's standard error
# and interval, vcov() of the log sd.
#   Rscript dev/quadseam-boundary.R <lib|base>
# log dev/quadseam-log/boundary-<tag>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
set.seed(3)
n <- 120
d1 <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(n))
d1$y <- rbinom(n, 1, plogis(0.2 + 0.5 * d1$x))
for (quad in c(FALSE, TRUE)) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(
    frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1,
        quadrature = quad),
    warning = function(x) {
      w <<- c(w, paste("WARN", class(x)[1], conditionMessage(x)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      w <<- c(w, paste("MSG", class(m)[1], conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  secs <- proc.time()[["elapsed"]] - t0
  vc <- suppressWarnings(VarCorr(f))
  s <- vc[["g"]]$sd
  sdr <- suppressWarnings(frmtmb:::sdr_of(f))
  V <- sdr$cov.fixed
  cat(sprintf("CASE quadrature=%-5s logLik %.9f theta %.4f secs %.2f\n",
              quad, as.numeric(logLik(f)), f$estimates$theta[1], secs))
  for (x in w) cat("  ", substr(gsub("[[:space:]]+", " ", x), 1, 200), "\n")
  cat(sprintf("  sd %.4g Est.Error %.4g Q2.5 %.4g Q97.5 %.4g\n",
              s[1, "Estimate"], s[1, "Est.Error"], s[1, "Q2.5"],
              s[1, "Q97.5"]))
  cat(sprintf("  cov.fixed diag %s  se_lost %s\n",
              paste(format(diag(V), digits = 4), collapse = " "),
              paste(names(sdr$se_lost), sdr$se_lost,
                    sep = "=", collapse = " ")))
  cat(sprintf("  vcov(beta) diag %s\n",
              paste(format(diag(suppressWarnings(vcov(f))), digits = 6),
                    collapse = " ")))
  ci <- suppressWarnings(confint(f))
  cat("  confint rows:", paste(rownames(ci), collapse = " "), "\n")
  print(ci)
}
cat("END\n")
