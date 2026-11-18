# Lane quadesc: RTMB::pnorm(log.p = TRUE) value and derivative on the tape
# far in the lower tail, where the skew normal's log density reads it.
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("RTMB", format(packageVersion("RTMB")), "\n")
f <- function(x) RTMB::pnorm(x, log.p = TRUE)
tp <- RTMB::MakeTape(f, 0)
J <- tp$jacfun()
for (x in c(-5, -20, -37, -38, -38.5, -39, -40, -100, -1e3)) {
  cat(sprintf("x %8g  value %.6g  exact %.6g  deriv %.6g  exact %.6g\n", x,
              tp(x), pnorm(x, log.p = TRUE), J(x),
              exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))))
}
