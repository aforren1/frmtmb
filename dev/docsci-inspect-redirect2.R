# build_redirect is where the target URL is formed. Read it, and
# copy_assets, so the two candidate routes for the 104 stubs can be
# compared on fact rather than on recollection.
for (nm in c("build_redirect", "copy_assets")) {
  cat("=====", nm, "\n")
  print(get(nm, envir = asNamespace("pkgdown")))
  cat("\n")
}
