# Release review of 0.74.0, item 3: render one vignette on rellib-r12,
# output and intermediates under dev/review-rel074-log/render/.
#   Rscript dev/review-rel074-render.R <path to .Rmd> [libs;..]
a <- commandArgs(trailingOnly = TRUE)
libs <- if (length(a) >= 2L) strsplit(a[2], ";", fixed = TRUE)[[1]] else
  "C:/Users/adf44/source/r/rellib-r12"
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
src <- normalizePath(a[1], "/")
tag <- gsub("[/:]", "_", sub(paste0(root, "/"), "", src, fixed = TRUE))
if (length(a) >= 2L) tag <- paste0(tag, "-altlib")
out <- file.path(root, "dev/review-rel074-log/render", sub("[.]Rmd$", "", tag))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# render a copy, so that knitr writes nothing beside the source
file.copy(src, file.path(out, basename(src)), overwrite = TRUE)
for (f in list.files(dirname(src), full.names = TRUE)) {
  if (!grepl("[.]Rmd$", f) && !dir.exists(f)) file.copy(f, out, overwrite = TRUE)
}
t0 <- Sys.time()
r <- tryCatch({
  rmarkdown::render(file.path(out, basename(src)), output_dir = out,
                    intermediates_dir = out, quiet = TRUE,
                    envir = new.env(parent = globalenv()))
  "OK"
}, error = function(e) paste("ERROR", conditionMessage(e)))
w <- warnings()
cat("RENDER", tag, r, "| libs:", paste(.libPaths()[1], collapse = ";"),
    "| frmtmb", as.character(packageVersion("frmtmb")), "|",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
