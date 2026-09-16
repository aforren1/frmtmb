## Reviewer's private library. NO roxygenise: the worktree must not be
## modified by the review, and man/ and NAMESPACE are already current in
## the lane's tree.
## Usage: Rscript dev/asrev-install.R
LIB <- "C:/Users/adf44/source/r/asrev-lib"
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-argspell"
dirs <- c(wt,
          file.path(wt, "extensions",
                    c("frmtmb.sample", "frmtmb.coupling", "frmtmb.spline",
                      "frmtmb.latent", "frmtmb.eam", "frmtmb.learn",
                      "frmtmb.ode")))
for (d in dirs) {
  cat("\n==== install", basename(d), "====\n")
  r <- system2(file.path(R.home("bin"), "R"),
               c("CMD", "INSTALL", paste0("--library=", shQuote(LIB)),
                 "--no-multiarch", shQuote(normalizePath(d))),
               stdout = TRUE, stderr = TRUE)
  cat(tail(r, 6), sep = "\n")
  if (!is.null(attr(r, "status")) && attr(r, "status") != 0) {
    cat(r, sep = "\n")
    stop("install failed: ", d)
  }
}
cat("\nlibrary now holds:\n")
print(utils::installed.packages(lib.loc = LIB)[, "Version", drop = FALSE])
