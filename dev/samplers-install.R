# Install the two outside samplers into the lane library only.
#   Rscript dev/samplers-install.R
# stanli's runtime goes under the lane library through R_USER_CACHE_DIR,
# so nothing lands in the user's own cache.
lib <- "C:/Users/adf44/source/r/wt-samplers-lib"
Sys.setenv(R_USER_CACHE_DIR = file.path(lib, "_cache"))
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
have <- function(p) nzchar(system.file(package = p))
cran <- "https://cloud.r-project.org"

if (!have("stanli")) {
  install.packages("stanli", repos = "https://seantalts.r-universe.dev",
                   lib = lib, dependencies = FALSE, type = "binary")
}
library(stanli)
stanli_install()
cat("stanli runtime:", stanli_runtime_path(), stanli_available(), "\n")

# SparseNUTS imports that the user library lacks
need <- c("ellipse", "microbenchmark", "sparseinv")
need <- need[!vapply(need, have, NA)]
if (length(need)) {
  install.packages(need, repos = cran, lib = lib, dependencies = FALSE,
                   type = "binary")
}
if (!have("StanEstimators")) {
  install.packages("StanEstimators",
                   repos = c("https://andrjohns.r-universe.dev", cran),
                   lib = lib, dependencies = FALSE, type = "binary")
}
if (!have("SparseNUTS")) {
  remotes::install_github("noaa-afsc/SparseNUTS",
                          ref = "2f3f1626219afce68fa2da0d884d4f2dca138117",
                          lib = lib, dependencies = FALSE, upgrade = "never",
                          build_vignettes = FALSE)
}
for (p in c("stanli", "StanEstimators", "SparseNUTS", "ellipse",
            "microbenchmark", "sparseinv")) {
  d <- packageDescription(p, lib.loc = lib)
  cat(p, d$Version, d$RemoteSha %||% d$Packaged %||% "", d$Repository %||% "",
      "\n")
}
