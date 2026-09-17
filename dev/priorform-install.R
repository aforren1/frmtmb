# Roxygenise the core package and install it into the lane's private
# library, and nothing else. The shared libraries are read-only.
#   Rscript dev/priorform-install.R [pkgdir]
lib <- "C:/Users/adf44/source/r/priorform-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
args <- commandArgs(trailingOnly = TRUE)
pkg <- if (length(args)) args[1] else "."
roxygen2::roxygenise(pkg)
cmd <- paste0("\"", file.path(R.home("bin"), "R.exe"), "\" CMD INSTALL ",
              "--no-multiarch --library=\"", lib, "\" \"", pkg, "\"")
status <- system(cmd)
if (status != 0) stop("install failed with status ", status)
cat("installed", pkg, "into", lib, "\n")
