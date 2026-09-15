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

# refit: importing generics' generic does NOT reach lme4's refit.merMod,
# because lme4 DEFINES its own refit rather than importing it, so its
# method lives in lme4's table.  Measured: dev/generics-out/fix-L.txt
# still lost refit after the import.  lme4 therefore joins the run-time
# table, with generics' generic as the fallback when lme4 is absent.
sub1("R/generic-owners.R",
  "  ngrps = c(\"brms\", \"lme4\")\n)",
  "  ngrps = c(\"brms\", \"lme4\"),\n  refit = \"lme4\"\n)")

sub1("R/generic-owners.R",
paste0("frm_adopt_generics <- function(pkgname = \"frmtmb\") {\n",
       "  owners <- unique(unlist(frm_generic_owners, use.names = FALSE))\n",
       "  adopted <- character()\n"),
paste0("frm_adopt_generics <- function(pkgname = \"frmtmb\") {\n",
       "  owners <- unique(unlist(frm_generic_owners, use.names = FALSE))\n",
       "  ns <- asNamespace(pkgname)\n",
       "  # `refit` reaches this namespace through importFrom(), so it\n",
       "  # has no binding of its own to switch later.  The namespace is\n",
       "  # still unsealed here and will not be again, so the shadow has\n",
       "  # to be made now, whether or not an owner ever turns up.\n",
       "  for (gen in names(frm_generic_owners)) {\n",
       "    if (exists(gen, envir = ns, inherits = FALSE)) next\n",
       "    g <- tryCatch(get(gen, envir = ns), error = function(e) NULL)\n",
       "    if (is.function(g)) assign(gen, g, envir = ns)\n",
       "  }\n",
       "  adopted <- character()\n"))

sub1("R/sugar.R",
paste0("#' @rdname refit\n#' @export\n",
       "refit.frmtmb_fit <- function(object, newresp, start = NULL, ...) {"),
paste0("#' @rdname refit\n#' @exportS3Method lme4::refit\n#' @export\n",
       "refit.frmtmb_fit <- function(object, newresp, start = NULL, ...) {"))
cat("DONE\n")
