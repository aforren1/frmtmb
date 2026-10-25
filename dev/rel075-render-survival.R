# Render frmtmb.survival's vignette from a copy outside the tree, on
# rellib-r13, and print the bhaz() section's output; then
# pkgdown::check_pkgdown() on frmtmb.survival (release fix 4).
#   Rscript dev/rel075-render-survival.R <scratch dir>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(RSTUDIO_PANDOC =
             "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
src <- "C:/Users/adf44/source/r/frmtmb-wt-rel075/extensions/frmtmb.survival"
dir.create(a[1], showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(src, "vignettes", "royston-parmar.Rmd"), a[1],
          overwrite = TRUE)
out <- rmarkdown::render(file.path(a[1], "royston-parmar.Rmd"),
                         quiet = TRUE)
cat("rendered:", out, file.size(out), "bytes\n")
md <- rmarkdown::render(file.path(a[1], "royston-parmar.Rmd"),
                        output_format = rmarkdown::md_document(),
                        quiet = TRUE)
x <- readLines(md)
i <- grep("spelling of a Cox baseline", x, fixed = TRUE)[1]
cat(x[i:length(x)], sep = "\n")
r <- tryCatch({ pkgdown::check_pkgdown(src); "OK" },
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("PKGDOWN frmtmb.survival", r, "\n")
