# Reviewer scratch: does R CMD check --as-cran set _R_CHECK_SUGGESTS_ONLY_,
# and where is the setting read?
x <- deparse(tools:::.check_packages)
j <- grep("as_cran", x, fixed = TRUE)
for (k in j) {
  blk <- x[k:min(length(x), k + 40)]
  hit <- grep("SUGGESTS_ONLY|DEPENDS_ONLY", blk)
  if (length(hit)) cat("as_cran block at", k, ":\n", blk[hit], sep = "\n")
}
cat("---- every SUGGESTS_ONLY / suggests_only line\n")
i <- grep("SUGGESTS_ONLY|suggests_only", x)
for (k in i) cat(k, ":", x[k], "\n")
