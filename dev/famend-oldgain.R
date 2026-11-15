# Lane famend, fix round F2: the new state assertions of
# test-family-fitend.R seen failing with the first round's
# frechet_ridge_gain(), which evaluated the probe and put the objective
# back without restoring TMB's best point. Swapped into this process only.
#   Rscript dev/famend-oldgain.R
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
g <- get("frechet_ridge_gain", ns)
b <- deparse(body(g))
i <- grep("saved <- obj_state_save", b)
b <- b[-c(i, i + 1L)]
b <- append(b, "  fit$obj$fn(p)", after = grep("f1 <- fit\\$obj\\$fn\\(q\\)", b))
body(g) <- parse(text = paste(b, collapse = "\n"))[[1]]
assignInNamespace("frechet_ridge_gain", g, ns)
print(body(get("frechet_ridge_gain", ns))[-(2:14)])
res <- test_file("tests/testthat/test-family-fitend.R", package = "frmtmb",
                 env = test_env("frmtmb"), reporter = "silent")
df <- as.data.frame(res)
for (r in res) for (e in r$results) {
  if (inherits(e, c("expectation_failure", "expectation_error"))) {
    cat("FAIL:", substr(gsub("\\s+", " ", conditionMessage(e)), 1, 200), "\n")
  }
}
cat(sprintf("RESULT pass=%d fail=%d error=%d\n", sum(df$passed),
            sum(df$failed), sum(df$error)))
