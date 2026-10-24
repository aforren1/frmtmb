# Render the latent vignette on the lane's build, into a scratch dir.
#   Rscript dev/mvnlatent-render.R <out dir>
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
o <- commandArgs(trailingOnly = TRUE)[1]
dir.create(o, showWarnings = FALSE)
file.copy("extensions/frmtmb.latent/vignettes/latent.Rmd", o,
          overwrite = TRUE)
rmarkdown::render(file.path(o, "latent.Rmd"), quiet = TRUE)
k <- knitr::purl(file.path(o, "latent.Rmd"), file.path(o, "latent.R"),
                 quiet = TRUE)
cat("rendered and tangled\n")
