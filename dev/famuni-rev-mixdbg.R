# Reviewer debug: does a NaN component density reach a mixture's row?
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("logspace_add(0, NaN) numeric:",
    tryCatch(RTMB::logspace_add(0, NaN), error = function(e) "err"), "\n")
f <- function(p) RTMB::logspace_add(p[1], p[2] + 0 * log(p[3]))
tp <- RTMB::MakeTape(f, c(0, 1, 1))
cat("on the tape, logspace_add(0, 1 + 0 * log(-1)):", tp(c(0, 1, -1)), "\n")
cat("jacobian:", tp$jacobian(c(0, 1, -1)), "\n")
f2 <- function(p) RTMB::logspace_add(p[1], p[2])
tp2 <- RTMB::MakeTape(f2, c(0, 1))
cat("on the tape, logspace_add(0, NaN):", tp2(c(0, NaN)),
    " logspace_add(NaN, 0):", tp2(c(NaN, 0)), "\n")
cat("jacobian at (0, NaN):", tp2$jacobian(c(0, NaN)), "\n")
