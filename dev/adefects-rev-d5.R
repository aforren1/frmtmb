source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb.latent))
cat("frmtmb.latent", format(packageVersion("frmtmb.latent")),
    "from", dirname(system.file(package = "frmtmb.latent")), "\n")
fns <- grep("^(hmm|lca)_", ls(asNamespace("frmtmb.latent")), value = TRUE)
fns <- fns[vapply(fns, function(f)
  "frmtmb.latent" %in% loadedNamespaces() &&
    f %in% getNamespaceExports("frmtmb.latent"), NA)]
cat("exported hmm_/lca_ readers:", paste(fns, collapse = ", "), "\n")
bad <- 0L; tot <- 0L
for (f in fns) {
  for (arg in list(1, "a", NULL, list())) {
    tot <- tot + 1L
    r <- tryCatch(do.call(f, list(arg)), error = function(e) e)
    if (!inherits(r, "frmtmb_error")) {
      bad <- bad + 1L
      cat("  UNCLASSED", f, "(", class(arg), "):",
          if (inherits(r, "condition"))
            substr(conditionMessage(r), 1, 70) else "NO ERROR", "\n")
    }
  }
}
cat(sprintf("classed frmtmb_error: %d of %d, unclassed %d\n",
            tot - bad, tot, bad))
r <- tryCatch(frmtmb.latent::hmm_starts(1), error = function(e) e)
cat("hmm_starts(1) class:", paste(class(r), collapse = ","), "\n")
cat("hmm_starts(1) msg:  ", substr(conditionMessage(r), 1, 110), "\n")
