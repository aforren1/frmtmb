# Which mechanism does the frmtmb in a given library carry?  Used to
# confirm that the reviewer's library still holds the SWAP build, so
# that the swap and the active binding can be timed in ONE interleaved
# run rather than compared across runs.
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib                     ", find.package("frmtmb"), "\n")
cat("frm_adopt_generics (swap)", exists("frm_adopt_generics", envir = ns,
                                       inherits = FALSE), "\n")
cat("frm_install_generics (ab)", exists("frm_install_generics",
                                        envir = ns, inherits = FALSE),
    "\n")
cat("loo binding is active   ", bindingIsActive("loo", ns), "\n")
cat("as_draws.frmtmb_fit     ",
    !is.null(tryCatch(getS3method("as_draws", "frmtmb_fit",
                                  optional = TRUE),
                      error = function(e) NULL)), "\n")
