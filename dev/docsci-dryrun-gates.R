# The gate list dev/release/build-docs.R derives, printed with the
# installed state of each, so the findings can name them.
root <- "."
rmd <- c(Sys.glob(file.path(root, "vignettes", "*.Rmd")),
         Sys.glob(file.path(root, "extensions", "*", "vignettes", "*.Rmd")))
cat("article sources scanned: ", length(rmd), "\n", sep = "")
txt <- unlist(lapply(rmd, readLines, warn = FALSE))
hit <- unlist(regmatches(txt, gregexpr('requireNamespace[(]"[^"]+"', txt)))
g <- sort(unique(sub('"$', "", sub('^requireNamespace[(]"', "", hit))))
for (p in g) {
  cat(sprintf("  %-14s %s\n", p,
              if (requireNamespace(p, quietly = TRUE)) "installed" else
                "MISSING"))
}
cat("gates: ", length(g), "\n", sep = "")
