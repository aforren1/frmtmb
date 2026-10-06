# Reviewer: precedence. A family's fit-end warning about one parameter
# (thres(): a threshold above every row) silences the SE warning even
# when the lost SE belongs to a different parameter (x, held by a bound).
# Also check_olre = "ignore": the structural check is silenced by the
# user, yet it still marks the fit explained.
#   Rscript dev/nanse-rev-prec.R base|lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
show <- function(lab, expr) {
  w <- character()
  fit <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  wv <- character()
  V <- withCallingHandlers(vcov(fit, full = TRUE), warning = function(x) {
    wv <<- c(wv, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  se <- suppressWarnings(sqrt(diag(V)))
  cat("\n==", lab, "\n  code", fit$opt$convergence, "\n")
  cat("  SE:", paste(sprintf("%s=%.3g", names(se), se), collapse = " "),
      "\n")
  cat("  lost:", paste(sprintf("%s(%s)", names(ns$sdr_of(fit)$se_lost),
                               ns$sdr_of(fit)$se_lost), collapse = " "),
      "\n")
  for (x in w) cat("  frm() warn:", substr(x, 1, 160), "\n")
  for (x in wv) cat("  vcov() warn:", substr(x, 1, 160), "\n")
  out <- capture.output(print(summary(fit)))
  i <- grep("without a standard error", out)
  if (length(i)) cat(paste0("  summary| ", out[i:min(length(out), i + 2)]),
                     sep = "\n")
}
set.seed(9)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
lat <- 1.0 * d$x + 0.5 * d$z + rlogis(n)
d$y <- 1L + (lat > -1) + (lat > 0) + (lat > 1)
show("cumulative, thres(4) on 4 categories, b held by ub = 0.1",
     frm(bf(y | thres(4) ~ x + z), family = cumulative(), data = d,
         prior = set_prior("", class = "b", ub = 0.1)))
show("control: same bound, thres() matching the data",
     frm(bf(y ~ x + z), family = cumulative(), data = d,
         prior = set_prior("", class = "b", ub = 0.1)))
set.seed(101)
d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
show("gaussian, (1 | g1) with one level, x held by ub = 0.1",
     frm(bf(y ~ x + (1 | g1)), family = gaussian(), data = d3,
         prior = set_prior("", class = "b", ub = 0.1)))
show("control: the same bound without (1 | g1)",
     frm(bf(y ~ x), family = gaussian(), data = d3,
         prior = set_prior("", class = "b", ub = 0.1)))
set.seed(5)
d2 <- data.frame(id = factor(1:80), x = rnorm(80))
d2$y <- 1 + 0.5 * d2$x + rnorm(80)
show("OLRE with check_olre = 'ignore'",
     frm(bf(y ~ x + (1 | id)), family = gaussian(), data = d2,
         control = frmtmb_control(check_olre = "ignore")))
