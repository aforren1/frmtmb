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
p <- "extensions/frmtmb.sample/R/methods-draws.R"

# R CMD check on frmtmb.sample: "Undocumented arguments in Rd file
# 'draws-dimensions.Rd': '...'", because nvariables() now carries the
# dots that posterior's generic has.
sub1(p,
paste0("#' @param x A `frmtmb_draws` from [frm_sample()].\n",
       "#' @return A single integer.\n"),
paste0("#' @param x A `frmtmb_draws` from [frm_sample()].\n",
       "#' @param ... Unused. `nvariables()` carries it because\n",
       "#'   posterior's generic does.\n",
       "#' @return A single integer.\n"))

sub1(p,
paste0("# posterior's nchains()/ndraws()/niterations()/nvariables() generics take\n",
       "# x alone, so these methods do too"),
paste0("# posterior's nchains(), ndraws() and niterations() generics take x\n",
       "# alone and nvariables() takes (x, ...); dev/generics-audit2.R in\n",
       "# core measured it. The signatures have to agree, because core now\n",
       "# hands these generics back to posterior."))
cat("DONE\n")
