# Lane optrap: does the standard-error check call a gaussian fit's
# coefficients flat when the response is in units 1e6 times larger?
# y ~ x (no random effects, which every build fits to the maximum) and
# the ri20 design, seeds 1 to 10.
#   Rscript dev/optrap-se1e6.R <lib or "">
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
for (form in c("y ~ x", "y ~ x + (1 | g)")) for (k in c(1, 1e3, 1e6)) {
  nw <- 0L
  short <- 0L
  for (s in 1:10) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- (1 + 0.5 * d$x + rnorm(100)) * k
    w <- character()
    f <- withCallingHandlers(frm(as.formula(form), data = d),
      warning = function(x) {w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")},
      message = function(x) invokeRestart("muffleMessage"))
    nw <- nw + any(grepl("Standard errors are not available", w))
    if (form == "y ~ x") {
      ll <- as.numeric(logLik(lm(y ~ x, data = d)))
      short <- short + (ll - as.numeric(logLik(f)) > 1e-6)
      gap <- max(c(if (s > 1) gap, ll - as.numeric(logLik(f))))
    }
  }
  cat(sprintf("%-16s x %-6g: SE warning on %d of 10%s\n", form, k, nw,
              if (form == "y ~ x") sprintf(", short of lm() on %d, max %.3g", short, gap)
              else ""))
}
