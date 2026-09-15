root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl(old, txt, fixed = TRUE)) stop("no match")
  if (length(gregexpr(old, txt, fixed = TRUE)[[1]]) != 1L) stop("many")
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

sub1("R/generic-owners.R",
paste0("#' Generics frmtmb re-exports from the package that owns them.\n",
       "#'\n",
       "#' nlme is Recommended and generics depends on base packages alone, so\n",
       "#' these four cost nothing and need no run-time machinery.\n"),
paste0("#' Generics frmtmb re-exports from the package that owns them.\n",
       "#'\n",
       "#' nlme is Recommended and ships with R; generics declares only base\n",
       "#' packages, `library(generics)` costing 0.000 s at the minimum of 12\n",
       "#' replicates against 0.860 s for lme4 (`dev/generics-loadcost.R`).\n",
       "#' So these four are free. Three of them need nothing further;\n",
       "#' `refit` also joins the run-time table, because lme4 defines a\n",
       "#' rival `refit` that this import cannot reach.\n"))

sub1("R/generic-owners.R",
paste0("#' Assign over a binding that may be locked.\n",
       "#'\n",
       "#' A namespace is sealed after `.onLoad` returns, so the load hook path\n",
       "#' has to unlock first. The binding is relocked afterwards, and only a\n",
       "#' name that is already there is written, so this never invents an\n",
       "#' export.\n"),
paste0("#' Assign over a binding that may be locked.\n",
       "#'\n",
       "#' A namespace is sealed after `.onLoad` returns, so the load hook\n",
       "#' path has to unlock first. The binding is relocked afterwards, and\n",
       "#' only a name that is ALREADY there is written, so this can never\n",
       "#' invent an export or reach a package that is not frmtmb.\n",
       "#'\n",
       "#' `unlockBinding()` is on R CMD check's list of possibly unsafe\n",
       "#' calls and raises a NOTE. It is here for ONE case and no other: an\n",
       "#' owner whose namespace is loaded after frmtmb has attached, and\n",
       "#' which is never attached itself, so that the search path cannot\n",
       "#' put its generic above frmtmb's. On the `.onLoad` path the\n",
       "#' namespace is not sealed yet and nothing is unlocked.\n"))
cat("DONE\n")
