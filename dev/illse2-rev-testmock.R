# Re-check of lane illse2, finding 4: run test-se-reparam.R with the
# asymmetry guard off (se_reparam_asym = Inf) or with the old scheme
# (Richardson only on single-parameter columns), and list failures.
#   Rscript dev/illse2-rev-testmock.R <noguard|schemeE>
mode <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-illse2-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({ library(testthat); library(frmtmb) })
ns <- asNamespace("frmtmb")
if (mode == "noguard") assignInNamespace("se_reparam_asym", Inf, "frmtmb")
if (mode == "schemeE") {
  rh <- ns$reparam_hessian
  wrap <- function(obj, p, M, step = 1e-3, rich_all = FALSE) {
    rh(obj, p, structure(M, lin_pos = seq_len(nrow(M))), step = step,
       rich_all = rich_all)
  }
  assignInNamespace("reparam_hessian", wrap, "frmtmb")
}
res <- test_file("C:/Users/adf44/source/r/frmtmb-wt-illse2/tests/testthat/test-se-reparam.R",
                 package = "frmtmb", env = test_env("frmtmb"),
                 reporter = "silent")
df <- as.data.frame(res)
cat(sprintf("mode %s: pass %d fail %d error %d\n", mode, sum(df$passed),
            sum(df$failed), sum(df$error)))
print(df[df$failed > 0 | df$error, c("test", "failed", "error")])
for (t in res) for (x in t$results) {
  if (inherits(x, "expectation_failure") || inherits(x, "expectation_error"))
    cat("  ", t$test, ":", substr(conditionMessage(x), 1, 200), "\n")
}
