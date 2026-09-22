# Reviewer, priority 3: the vignettes. The lane found three failing
# only at R CMD build time, so a grep is not the instrument. This
# KNITS every vignette of all eight packages against the lane library,
# one R process per vignette, and reports the first error of each.
#
#   Rscript dev/shapes-rev-vig.R <path/to.Rmd>

arg <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
f <- normalizePath(arg[1], winslash = "/")
out <- tempfile(fileext = ".md")
t0 <- proc.time()
r <- tryCatch({
  # knit in the vignette's own directory, as R CMD build does
  old <- setwd(dirname(f)); on.exit(setwd(old))
  suppressWarnings(suppressMessages(
    knitr::knit(basename(f), output = out, quiet = TRUE)))
  "OK"
}, error = function(c) paste0("ERROR: ", conditionMessage(c)))
cat(sprintf("%-58s %7.1fs  %s\n", basename(f), (proc.time() - t0)[3],
            substr(gsub("[\r\n]+", " ", r), 1, 200)))
