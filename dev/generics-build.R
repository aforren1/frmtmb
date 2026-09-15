# Roxygenise and install core into the lane's private library.
LIB <- "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
args <- commandArgs(trailingOnly = TRUE)
pkg <- if (length(args)) args[1] else root

cat("roxygen2", as.character(packageVersion("roxygen2")), "\n")
roxygen2::roxygenise(pkg, load_code = roxygen2::load_pkgload)
cat("--- roxygenised ---\n")

o <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", paste0("--library=", shQuote(LIB)),
               "--no-multiarch", shQuote(pkg)),
             stdout = TRUE, stderr = TRUE)
cat(paste(o, collapse = "\n"), "\n")
cat("--- installed ---\n")
