m <- parseNamespaceFile("frmtmb-wt-generics", "C:/Users/adf44/source/r")$S3methods
print(dim(m))
print(m[m[, 1] %in% c("as_draws", "ngrps", "refit", "loo_compare", "posterior_summary"), ])
