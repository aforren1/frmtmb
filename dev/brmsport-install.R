# Install this worktree's frmtmb and frmtmb.sample into the lane's private
# library ONLY. The shared reference build rellib-r3 was built at 02:17
# UTC on 2026-09-17, before commit 4e179c0 last touched R/ at 02:58 UTC,
# so it is not provably the base commit.
#
#   Rscript dev/brmsport-install.R [frmtmb] [frmtmb.sample]
lib <- "C:/Users/adf44/source/r/brmsport-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(identical(normalizePath(.libPaths()[1], winslash = "/"),
                    normalizePath(lib, winslash = "/")))
root <- "C:/Users/adf44/source/r/frmtmb-wt-brmsport"
pkgs <- commandArgs(trailingOnly = TRUE)
if (!length(pkgs)) pkgs <- c("frmtmb", "frmtmb.sample")
for (p in pkgs) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  st <- system2(file.path(R.home("bin"), "R.exe"),
                c("CMD", "INSTALL", "--no-multiarch",
                  paste0("--library=", shQuote(lib)), shQuote(path)))
  if (st != 0) stop("install of ", p, " failed")
}
cat("installed", pkgs, "into", lib, "\n")
