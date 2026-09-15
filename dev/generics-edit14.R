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

# Found by its own benchmark: calling frm_adopt_generics() twice
# APPENDED a second copy of every load hook, because setHook() appends.
# .onLoad runs once per session so nothing shipped was wrong, but
# devtools::load_all() and any re-load would stack them, and the
# benchmark that measured the loop was measuring a growing hook list.
# Dropping the previous hooks first makes the function idempotent.
sub1("R/generic-owners.R",
paste0("frm_adopt_generics <- function(pkgname = \"frmtmb\") {\n",
       "  owners <- unique(unlist(frm_generic_owners, use.names = FALSE))\n",
       "  ns <- asNamespace(pkgname)\n"),
paste0("frm_adopt_generics <- function(pkgname = \"frmtmb\") {\n",
       "  owners <- unique(unlist(frm_generic_owners, use.names = FALSE))\n",
       "  ns <- asNamespace(pkgname)\n",
       "  # setHook() APPENDS, so a second call would register a second\n",
       "  # copy of every hook. .onLoad runs once per session, but\n",
       "  # devtools::load_all() and a re-load do not, and the benchmark\n",
       "  # that timed this loop was measuring a growing hook list.\n",
       "  frm_drop_generic_hooks()\n"))
cat("DONE\n")
