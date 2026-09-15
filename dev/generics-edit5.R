root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (n != 1L) stop(n, " matches in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

sub1("R/generic-owners.R",
paste0("#  * Static. `fixef`, `ranef` and `VarCorr` come from nlme, which is a\n",
       "#    Recommended package that ships with R, and `refit` from generics,\n",
       "#    which depends on base packages alone. Both are importFrom()-ed and\n",
       "#    re-exported, which is exactly what lme4 and brms do with nlme.\n",
       "#    Nothing is left to run time.\n"),
paste0("#  * Static. `fixef`, `ranef` and `VarCorr` come from nlme, which is a\n",
       "#    Recommended package that ships with R, and `refit` from generics,\n",
       "#    which depends on base packages alone. Both are importFrom()-ed and\n",
       "#    re-exported, which is exactly what lme4 and brms do with nlme.\n",
       "#    nlme is then the ONE generic that frmtmb, lme4, glmmTMB and brms\n",
       "#    all dispatch through, because all three of the others import it\n",
       "#    from nlme too.\n",
       "#\n",
       "#    `refit` is the exception in that group and needs both routes.\n",
       "#    lme4 DEFINES its own `refit` rather than importing generics', so\n",
       "#    `refit.merMod` lives in lme4's table and the generics import does\n",
       "#    not reach it: measured in dev/generics-out/fix-L.txt, `refit` was\n",
       "#    the one name still lost to an lme4 user after the import landed.\n",
       "#    So lme4 is in the run-time table for `refit` as well, with\n",
       "#    generics' generic as the fallback when lme4 is not there.\n"))
cat("DONE\n")
