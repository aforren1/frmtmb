p <- "R/autocor.R"
txt <- paste(readLines(p), collapse = "\n")
old <- "  if (autocor_is_seam(ac)) return(ac_block_entry(ac)[[\"start\"]](ac[[\"aux\"]]))\n"
new <- "  if (autocor_is_seam(ac)) {\n    return(ac_block_entry(ac)[[\"start\"]](ac[[\"aux\"]]))\n  }\n"
stopifnot(grepl(old, txt, fixed = TRUE))
txt <- sub(old, new, txt, fixed = TRUE)
old <- "  if (autocor_is_seam(ac)) return(rep(\"raw\", length(autocor_natural(\n    numeric(ac[[\"npar\"]]), ac))))\n"
new <- "  if (autocor_is_seam(ac)) {\n    return(rep(\"raw\", length(autocor_natural(numeric(ac[[\"npar\"]]), ac))))\n  }\n"
stopifnot(grepl(old, txt, fixed = TRUE))
txt <- sub(old, new, txt, fixed = TRUE)
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
