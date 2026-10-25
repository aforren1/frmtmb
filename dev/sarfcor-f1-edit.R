# Fix round 1, part A: the seam's bounds slot. Rscript <this> <root>
root <- commandArgs(TRUE)[1]
p <- file.path(root, "R/autocor-seam.R")
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("                                    natural = NULL, sigma_varies = TRUE,
                                    usage = NULL, compat_note = NULL,",
"                                    natural = NULL, sigma_varies = TRUE,
                                    bounds = NULL, usage = NULL,
                                    compat_note = NULL,")
rep1("  natural <- fn_or(natural, \"natural\", 2L, null_ok = TRUE)\n",
"  natural <- fn_or(natural, \"natural\", 2L, null_ok = TRUE)
  bounds <- fn_or(bounds, \"bounds\", 1L, null_ok = TRUE)\n")
rep1("                natural = natural, sigma_varies = sigma_varies,
                usage = usage, package = package)",
"                natural = natural, sigma_varies = sigma_varies,
                bounds = bounds, usage = usage, package = package)")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
