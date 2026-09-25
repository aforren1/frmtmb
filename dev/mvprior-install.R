# Lane wt-mvprior: roxygenise and install packages into the lane library
# only. Usage: Rscript dev/mvprior-install.R core [sample] [<ext> ...]
# where <ext> is an extension directory name such as frmtmb.eam.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
which <- commandArgs(trailingOnly = TRUE)
for (w in which) {
  d <- if (identical(w, "core")) MVPRIOR_ROOT else
    file.path(MVPRIOR_ROOT, "extensions",
              if (identical(w, "sample")) "frmtmb.sample" else w)
  cat("== roxygenise", d, "\n")
  roxygen2::roxygenise(d)
  cat("== install", d, "\n")
  st <- system2(file.path(R.home("bin"), "R.exe"),
                c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
                  paste0("--library=", MVPRIOR_LIB), shQuote(d)))
  if (st != 0) stop("install failed: ", d)
}
print(installed.packages(lib.loc = MVPRIOR_LIB)[, "Version"])
