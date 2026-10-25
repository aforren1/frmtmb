# What R CMD check --as-cran sets that decides which installed packages
# a package's tests can see: read from this R's own tools sources.
#   Rscript dev/smallfix-ascran-env.R
# Log: dev/smallfix-out/ascran-env.txt
b <- deparse(tools:::.check_packages)
cat(R.version.string, "\n")
i <- grep("SUGGESTS_ONLY|INSTALL_DEPENDS|R_check_suggests_only|install_depends",
          b)
show <- sort(unique(unlist(lapply(i, function(k) max(1, k - 3):(k + 6)))))
for (k in show) cat(k, ":", trimws(b[k]), "\n")
