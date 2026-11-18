# Review of lane quadesc: RTMB pnorm(log.p = TRUE) value and derivative
# beyond the lane's x = -1000, where a run-off alpha (1e6 to 1e7) puts it.
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("RTMB", format(packageVersion("RTMB")), "\n")
f <- function(x) RTMB::pnorm(x, log.p = TRUE)
tp <- RTMB::MakeTape(f, 0)
J <- tp$jacfun()
for (x in c(-1e3, -1e4, -1e5, -1e6, -1e7, -1e8, -1e10, -1e100, -1e154,
            -1e155, 1e3, 1e6, 1e8)) {
  cat(sprintf("x %10g  value %.6g  exact %.6g  deriv %.6g\n", x,
              tp(x), pnorm(x, log.p = TRUE), J(x)))
}
