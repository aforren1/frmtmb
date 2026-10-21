# A quadrature fit with a variance at zero: what do the SE check,
# VarCorr() and vcov() report? One block (r10 and lane alike) and the
# nested fuzz spec (lane only; r10's summary() fails there).
#   Rscript dev/fuzzfix-quad-boundary.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(3)
n <- 120
d1 <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(n))
d1$y <- rbinom(n, 1, plogis(0.2 + 0.5 * d1$x))
d2 <- readRDS("dev/fuzzfix-out/data-cum.rds")
cases <- list(
  one = list(bf(y ~ x + (1 | g)) + bernoulli(), d1),
  one_lap = list(bf(y ~ x + (1 | g)) + bernoulli(), d1, FALSE),
  nested = list(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) + cumulative(),
                d2),
  nested_lap = list(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) +
                      cumulative(), d2, FALSE))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  quad <- if (length(cs) > 2) cs[[3]] else TRUE
  w <- character()
  f <- withCallingHandlers(frm(cs[[1]], data = cs[[2]], quadrature = quad),
    warning = function(x) { w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning") },
    message = function(m) { w <<- c(w, paste("MSG", conditionMessage(m)))
      invokeRestart("muffleMessage") })
  if (nm == "nested" && a[1] == "base") next
  vc <- suppressWarnings(VarCorr(f))
  th <- grep("^theta|^sd_", rownames(vcov(f)), value = TRUE)
  cat(sprintf("CASE %-10s warn/msg: %s\n", nm,
              paste(substr(w, 1, 70), collapse = " || ")))
  for (g in names(vc)) {
    s <- vc[[g]]$sd
    cat(sprintf("  %-6s sd %.4g Est.Error %.4g Q2.5 %.4g Q97.5 %.4g\n", g,
                s[1, "Estimate"], s[1, "Est.Error"], s[1, "Q2.5"],
                s[1, "Q97.5"]))
  }
}
