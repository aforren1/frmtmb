source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")

# N1. Item 2.5d's whole deliverable is a page the user is told to read,
# and it was the only one of the nine comparable frmtmb-* concept
# pages that was @keywords internal and absent from _pkgdown.yml.
sub1("R/scales.R",
     "#' @name frmtmb-scales\n#' @keywords internal\nNULL\n",
     "#' @name frmtmb-scales\nNULL\n")
sub1("_pkgdown.yml",
     "  - sigma.frmtmb_fit\n",
     "  - sigma.frmtmb_fit\n  - frmtmb-scales\n")
cat("DONE\n")
