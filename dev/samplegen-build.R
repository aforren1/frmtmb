# Roxygenise and install one package of this worktree into the lane's
# private library, ONLY that library.
#   Rscript dev/samplegen-build.R <package dir> [library]
av <- commandArgs(trailingOnly = TRUE)
LIB <- if (length(av) >= 2L) av[[2L]] else
  "C:/Users/adf44/source/r/samplegen-lib"
dir.create(LIB, showWarnings = FALSE)
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkg <- av[[1L]]
cat("roxygen2", as.character(packageVersion("roxygen2")), "\n")
roxygen2::roxygenise(pkg, load_code = roxygen2::load_pkgload)
cat("--- roxygenised ---\n")
o <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", paste0("--library=", shQuote(LIB)),
               "--no-multiarch", shQuote(pkg)),
             stdout = TRUE, stderr = TRUE)
cat(paste(o, collapse = "\n"), "\n")
cat("--- installed into", LIB, "---\n")
