## Evidence for the thirteen dots exemptions: the generics that emmeans,
## insight and marginaleffects own are called BY THOSE PACKAGES with
## arguments of the caller's choosing, so a method that refuses an
## unknown name there breaks the caller and not the caller's typo.
##
## Run: Rscript dev/argspell-exempt-evidence.R
.libPaths(c("C:/Users/adf44/source/r/argspell-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

show_calls <- function(pkg, pattern) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, ": NOT INSTALLED\n"); return(invisible(NULL))
  }
  ns <- asNamespace(pkg)
  hits <- 0L
  for (nm in ls(ns, all.names = TRUE)) {
    o <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
    if (!is.function(o)) next
    src <- tryCatch(deparse(body(o)), error = function(e) character())
    ln <- grep(pattern, src, value = TRUE)
    ln <- grep("[.][.][.]", ln, value = TRUE)
    if (length(ln)) {
      hits <- hits + length(ln)
      cat(sprintf("  %s::%s  %s\n", pkg, nm, trimws(ln[1L])))
    }
  }
  cat(sprintf("%s: %d call sites forwarding dots into %s\n\n",
              pkg, hits, pattern))
}

show_calls("emmeans", "emm_basis\\(")
show_calls("emmeans", "recover_data\\(")
show_calls("insight", "find_formula\\(")
show_calls("insight", "get_parameters\\(")
show_calls("marginaleffects", "get_predict\\(")
