# Reviewer mutation helper. mut() rewrites one function of a namespace
# by a literal pattern on its deparsed text and refuses when the
# pattern does not match exactly `times` times, so a mutant that
# silently changes nothing cannot pass for a caught one.
esc <- function(s) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", s)
pat_of <- function(from) gsub("\\s+", "\\\\s*", esc(from))
mut <- function(fn, from, to, pkg = "frmtmb", times = 1L) {
  ns <- asNamespace(pkg)
  f <- get(fn, envir = ns)
  txt <- paste(deparse(f), collapse = "\n")
  p <- pat_of(from)
  hits <- gregexpr(p, txt, perl = TRUE)[[1]]
  k <- sum(hits > 0)
  if (k != times) stop("mutant ", fn, ": pattern matched ", k,
                       " times, expected ", times, ": ", from)
  txt2 <- gsub(p, to, txt, perl = TRUE)
  g <- eval(parse(text = txt2))
  environment(g) <- environment(f)
  put(fn, g, pkg)
}
put <- function(fn, g, pkg = "frmtmb") {
  ns <- asNamespace(pkg)
  utils::assignInNamespace(fn, g, ns = ns)
  pe <- paste0("package:", pkg)
  if (pe %in% search() && exists(fn, envir = as.environment(pe),
                                 inherits = FALSE)) {
    e <- as.environment(pe)
    unlockBinding(fn, e); assign(fn, g, envir = e); lockBinding(fn, e)
  }
  # every importer's copy, and the S3 table when fn is a method
  for (imp in loadedNamespaces()) {
    ie <- parent.env(asNamespace(imp))
    if (exists(fn, envir = ie, inherits = FALSE) &&
        identical(environment(get(fn, envir = ie)), ns)) {
      unlockBinding(fn, ie); assign(fn, g, envir = ie); lockBinding(fn, ie)
    }
  }
  if (grepl("[.]frmtmb_", fn)) {
    gen <- sub("[.]frmtmb_.*$", "", fn); cls <- sub("^.*?[.](frmtmb_)", "\\1",
                                                    fn, perl = TRUE)
    registerS3method(gen, cls, g, envir = ns)
    for (own in c("brms", "stats", "loo", "bayesplot", "posterior")) {
      if (isNamespaceLoaded(own) &&
          exists(gen, envir = asNamespace(own), inherits = FALSE)) {
        registerS3method(gen, cls, g, envir = asNamespace(own))
      }
    }
  }
  cat("MUTATED ", pkg, "::", fn, "\n", sep = "")
}
