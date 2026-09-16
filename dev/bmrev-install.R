## Reviewer's private install. Installs ONLY into bmrev-lib.
LIB <- "C:/Users/adf44/source/r/bmrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("libPaths:\n"); print(.libPaths())
pkg <- "C:/Users/adf44/source/r/frmtmb-wt-brmsmatch/extensions/frmtmb.sample"
cat("installing", pkg, "into", LIB, "\n")
res <- system2("R", c("CMD", "INSTALL", paste0("--library=", LIB), pkg),
               stdout = TRUE, stderr = TRUE)
cat(res, sep = "\n")
cat("\n== versions ==\n")
for (p in c("frmtmb", "frmtmb.sample", "StanHeaders", "rstan", "brms",
            "posterior", "bayesplot", "matrixStats")) {
  cat(sprintf("%-14s %-10s %s\n", p,
              tryCatch(format(packageVersion(p)), error = function(e) "?"),
              tryCatch(dirname(system.file(package = p)),
                       error = function(e) "?")))
}
cat("DONE\n")
