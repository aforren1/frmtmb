# build_sitemap is not exported. Read its signature and body, because
# the root sitemap also covers the subsites and a core-first build loses
# that the same way the search index does.
f <- get("build_sitemap", envir = asNamespace("pkgdown"))
print(args(f))
print(f)
