root <- commandArgs(trailingOnly = TRUE)[1]; p <- file.path(root, "R/par-template.R")
s <- readChar(p, file.info(p)$size, useBytes = TRUE)
a <- "if (from_flat && same(ll, floor)) break"; b <- "if (agrees && !(saddle_floor && same(ll, floor))) break"
stopifnot(grepl(a, s, fixed = TRUE), grepl(b, s, fixed = TRUE))
s <- sub(a, "if (Sys.getenv(\"REV_MUT\") != \"no_early\" && from_flat && same(ll, floor)) break", s, fixed = TRUE)
s <- sub(b, "if (Sys.getenv(\"REV_MUT\") != \"no_early\" && agrees && !(saddle_floor && same(ll, floor))) break", s, fixed = TRUE)
writeBin(charToRaw(s), p); cat("mutated\n")
