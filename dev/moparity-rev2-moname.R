# Reviewer, re-check: a data column NAMED mo (not a call to mo()) in an
# interaction and in a || term, base against lane.
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(2); n <- 300
d <- data.frame(mo = rnorm(n), z = rnorm(n), fm = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(sample(1:20, n, TRUE)))
d$mo2 <- d$fm
d$y <- 1 + 0.3 * d$mo * d$z + rnorm(20)[d$g] + rnorm(n)
for (fo in list(y ~ mo * z, y ~ mo:z, y ~ mo + (1 + mo || g), y ~ z + (1 + z || g))) {
  r <- tryCatch(suppressWarnings(frm(bf(fo), family = gaussian(), data = d)), error = function(e) e)
  cat(deparse1(fo), ":", if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else
    paste("logLik", format(as.numeric(logLik(r)), digits = 10), "| blocks",
          paste(vapply(r$frame$re_blocks, function(b) paste(unlist(lapply(b$components, `[[`, "cnms")), collapse = ","), ""), collapse = " ; ")), "\n")
}
