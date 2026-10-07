# Reviewer of the 0.69.0 consolidation: optima's c1 tests seen to fail
# without c1. summary_mo_frame()'s face rule is put back to a weight
# below mo_face_tol only (the punch-2 rule), in frmtmb's namespace and
# in frmtmb.sample's imports, and the test file is run.
#   Rscript dev/relrev069-c1-before.R <package> <test file>
.libPaths(c("C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
suppressMessages({library(testthat); library(a[1], character.only = TRUE)})
ns <- asNamespace("frmtmb")
f <- get("summary_mo_frame", ns)
src <- deparse(f)
i <- grep("face <- w < mo_face_tol | w > 1 - mo_face_tol", src, fixed = TRUE)
stopifnot(length(i) == 1L)
src[i] <- sub(" | w > 1 - mo_face_tol", "", src[i], fixed = TRUE)
g <- eval(parse(text = src)); environment(g) <- ns
unlockBinding("summary_mo_frame", ns); assign("summary_mo_frame", g, ns)
if (isNamespaceLoaded("frmtmb.sample")) {
  imp <- parent.env(asNamespace("frmtmb.sample"))
  if (exists("summary_mo_frame", imp, inherits = FALSE)) {
    unlockBinding("summary_mo_frame", imp)
    assign("summary_mo_frame", g, imp)
    cat("patched frmtmb.sample's import too\n")
  }
}
res <- as.data.frame(test_file(a[2], package = a[1],
                               env = test_env(a[1]), reporter = "silent"))
print(res[res$failed > 0 | res$error, c("test", "nb", "failed", "error")])
cat("RESULT pass=", sum(res$passed), " fail=", sum(res$failed),
    " err=", sum(res$error), "\n", sep = "")
