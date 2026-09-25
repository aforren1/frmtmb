.libPaths(c("C:/Users/adf44/source/r/predfix-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
tree <- "C:/Users/adf44/source/r/frmtmb-wt-predfix"
pk <- list(frmtmb = tree,
           frmtmb.sample = file.path(tree, "extensions/frmtmb.sample"),
           frmtmb.latent = file.path(tree, "extensions/frmtmb.latent"),
           frmtmb.coupling = file.path(tree, "extensions/frmtmb.coupling"),
           frmtmb.learn = file.path(tree, "extensions/frmtmb.learn"))
for (p in names(pk)) {
  ns <- asNamespace(p)
  cat(p, as.character(packageVersion(p)), find.package(p), "\n")
  env <- new.env(parent = ns)
  fs <- list.files(file.path(pk[[p]], "R"), "[.][Rr]$", full.names = TRUE)
  for (f in fs) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && (identical(e[[1]], as.name("<-")) ||
                         identical(e[[1]], as.name("="))) &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) {
        assign(as.character(e[[2]]), eval(e[[3]], env), env)
      }
    }
  }
  nms <- ls(env, all.names = TRUE)
  diffs <- character(0); miss <- character(0)
  for (n in nms) {
    if (!exists(n, envir = ns, inherits = FALSE)) { miss <- c(miss, n); next }
    a <- get(n, envir = ns, inherits = FALSE)
    if (!is.function(a)) next
    b <- get(n, envir = env)
    if (!identical(deparse(body(a)), deparse(body(b))) ||
        !identical(deparse(formals(a)), deparse(formals(b))))
      diffs <- c(diffs, n)
  }
  cat("  functions compared", length(nms), "differ", length(diffs),
      "missing", length(miss), "\n")
  if (length(diffs)) cat("  DIFF:", head(diffs, 30), "\n")
  if (length(miss)) cat("  MISSING:", head(miss, 30), "\n")
}
cat("cs_offsets_add exported:",
    "cs_offsets_add" %in% getNamespaceExports("frmtmb"), "\n")
