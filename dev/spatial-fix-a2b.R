# Fix round 1, A2: rewrap the OLRE comment fix-a2 wrote, on the tree
# named as the argument.
root <- commandArgs(trailingOnly = TRUE)[1]
p <- file.path(root, "R/fit.R")
s <- paste(readLines(p, warn = FALSE), collapse = "\n")
old <- "    # and the residual sd is not a reparameterization of it. That is the animal model: the additive\n    # genetic variance and the residual variance are separately\n    # identified precisely BECAUSE A is not the identity."
new <- "    # and the residual sd is not a reparameterization of it. That is the\n    # animal model: the additive genetic variance and the residual\n    # variance are separately identified precisely BECAUSE A is not the\n    # identity."
if (grepl(old, s, fixed = TRUE)) writeLines(sub(old, new, s, fixed = TRUE), p)
cat("fix-a2b applied to", root, "\n")
