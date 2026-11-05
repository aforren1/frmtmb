# Fix round 2's tests, each seen to fail with its fix taken out, in-process.
# Usage: Rscript dev/latres-mut.R <mutant>
m <- commandArgs(TRUE)[[1]]
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(testthat)
suppressPackageStartupMessages(library(frmtmb))
setv <- function(n, v) assignInNamespace(n, v, "frmtmb")
switch(m,
  none = NULL,
  r1 = {
    setv("latres_disp_upper", function(frame, spec, upper, map) upper)
    setv("latres_disp_report", function(fit) invisible(NULL))
  },
  margin = setv("latres_refuse_margin", function(fit) invisible(NULL)),
  r2_flat = setv("latres_se_flat", -1),
  r2_block = setv("latres_limit_par", function(fit, nms) grepl("^theta_", nms)),
  r3_cubic = setv("latres_cubic_floor", -1),
  r3_disp = {
    setv("latres_share_disp", 0)
    setv("latres_info_disp", 0)
  },
  r4 = {
    f <- frmtmb:::latres_skip_se
    setv("latres_skip_se", function(fit, newdata = NULL) f(fit, NULL))
  },
  r5 = setv("latres_gh_lost", 0),
  stop("no such mutant"))
res <- test_file("tests/testthat/test-latres.R", package = "frmtmb",
                 env = test_env("frmtmb"), reporter = "silent",
                 stop_on_failure = FALSE)
df <- as.data.frame(res)
bad <- df[df$failed > 0 | df$error, ]
cat(sprintf("MUT %s: tests=%d fail=%d err=%d\n", m, sum(df$nb), sum(df$failed), sum(df$error)))
for (i in seq_len(nrow(bad))) cat("  ", bad$test[i], "\n")
