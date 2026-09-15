## Targeted revalidation after punch round 1, in place of a second
## `R CMD check --as-cran`.
##
## What the punch round changed is documentation and one test: the Rd
## section, the vignette prose, NEWS, and an assertion that costs no
## fit. The check's expensive halves, examples and tests, are unchanged
## or rerun per file. So this runs the two parts that could actually
## break: the Rd has to parse and validate, and the vignette has to
## knit.
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
rd <- "extensions/frmtmb.coupling/man/cross_wishart.Rd"
p <- tools::parse_Rd(rd)
cat("parse_Rd OK, sections:", length(p), "\n")
out <- utils::capture.output(
  tools::checkRd(rd), type = "message")
cat("checkRd messages:", length(out), "\n")
if (length(out)) writeLines(out)
## Every Rd in the package, since roxygen rewrote from one source.
for (f in list.files("extensions/frmtmb.coupling/man", full.names = TRUE)) {
  o <- utils::capture.output(tools::checkRd(f), type = "message")
  if (length(o)) {
    cat("--", basename(f), "\n")
    writeLines(o)
  }
}
cat("all man pages checked\n")

## The vignette, knitted the way R CMD check rebuilds it.
Sys.setenv(RSTUDIO_PANDOC =
  "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
od <- file.path(tempdir(), "coh-vig")
dir.create(od, showWarnings = FALSE)
t0 <- Sys.time()
knitr::knit(input = "extensions/frmtmb.coupling/vignettes/coherence.Rmd",
            output = file.path(od, "coherence.md"), quiet = TRUE)
cat(sprintf("vignette knitted in %.0f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
