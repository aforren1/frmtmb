# Lane wt-mvprior: does the lane library hold this tree's code? For each
# package, every function defined in its R/ sources is compared with the
# installed namespace's, by deparsed body and formals.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
LIB <- Sys.getenv("LIBMATCH_LIB", MVPRIOR_LIB)
pk <- c(frmtmb = ".", frmtmb.sample = "extensions/frmtmb.sample",
        frmtmb.latent = "extensions/frmtmb.latent",
        frmtmb.eam = "extensions/frmtmb.eam")
for (p in names(pk)) {
  ns <- asNamespace(loadNamespace(p, lib.loc = LIB))
  stopifnot(startsWith(getNamespaceInfo(ns, "path"), LIB))
  same <- 0L
  diff <- character(0)
  for (f in list.files(file.path(MVPRIOR_ROOT, pk[[p]], "R"), "[.]R$",
                       full.names = TRUE)) {
    for (e in as.list(parse(f, keep.source = FALSE))) {
      if (!is.call(e) || !identical(e[[1]], as.name("<-")) ||
          !is.name(e[[2]]) || !is.call(e[[3]]) ||
          !identical(e[[3]][[1]], as.name("function"))) next
      nm <- as.character(e[[2]])
      if (!exists(nm, envir = ns, inherits = FALSE)) next
      inst <- get(nm, envir = ns)
      if (!is.function(inst)) next
      src <- eval(e[[3]], baseenv())
      ok <- identical(deparse(body(src)), deparse(body(inst))) &&
        identical(deparse(formals(src)), deparse(formals(inst)))
      if (ok) same <- same + 1L else diff <- c(diff, nm)
    }
  }
  cat(p, "identical", same, "differ", length(diff),
      if (length(diff)) paste0("(", paste(diff, collapse = ", "), ")"), "\n")
}
