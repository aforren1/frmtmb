# The release review's F2: under nl_starts = FALSE the saddle restart
# still runs and says so, so the text must not promise "no message about
# a restart". The scope stays the approved one (the sign probe and the
# pole restart); the text now says that the saddle restart still runs.
f <- "C:/Users/adf44/source/r/frmtmb-wt-rel083/R/fit.R"
x <- readLines(f)
i <- which(x == "#'   runs both. `FALSE` runs neither, so the fit is the one the optimizer")
stopifnot(length(i) == 1, x[i + 4L] == "#'   (`frmtmb_nl_limit`) is not affected.")
x <- c(x[seq_len(i - 1L)],
  "#'   runs both. `FALSE` turns off the sign probe and the pole restart. A",
  "#'   fit can still restart from a saddle point (message class",
  "#'   `frmtmb_nl_restart`); this switch does not change that, nor the",
  "#'   warning that a fit is on a ridge to a limit of its body",
  "#'   (`frmtmb_nl_limit`).",
  x[-seq_len(i + 4L)])
con <- file(f, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
cat("done\n")
