# Dry run of the "Install the core package's Suggests" step in
# .github/workflows/pkgdown.yaml. Same code, minus the install, against
# this machine's library so the skip logic is exercised too.
d <- read.dcf("DESCRIPTION")
p <- unlist(strsplit(d[, "Suggests"], ","))
p <- trimws(sub("[(].*", "", p))
p <- p[nzchar(p)]
p <- p[!startsWith(p, "frmtmb")]
cat("core Suggests parsed: ", length(p), "\n", sep = "")
miss <- setdiff(p, rownames(installed.packages()))
cat("would install here: ", length(miss), " -> ",
    paste(miss, collapse = " "), "\n", sep = "")
cat("all parsed names:\n  ", paste(p, collapse = " "), "\n", sep = "")
