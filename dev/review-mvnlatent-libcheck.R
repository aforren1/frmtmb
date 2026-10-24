# Is the lane library the worktree source? Deparse every function in
# the source R/ files and compare with the installed namespace.
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-mvnlatent"
dirs <- c(frmtmb = file.path(wt, "R"),
          frmtmb.latent = file.path(wt, "extensions/frmtmb.latent/R"),
          frmtmb.sample = file.path(wt, "extensions/frmtmb.sample/R"))
for (p in names(dirs)) {
  ns <- asNamespace(p)
  env <- new.env(parent = ns)
  for (f in sort(list.files(dirs[[p]], "\\.R$", full.names = TRUE))) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && (identical(e[[1]], as.name("<-")) ||
                         identical(e[[1]], as.name("="))) &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) {
        assign(as.character(e[[2]]), eval(e[[3]], env), envir = env)
      }
    }
  }
  nm <- ls(env, all.names = TRUE)
  dif <- character(0)
  miss <- character(0)
  for (n in nm) {
    if (!exists(n, envir = ns, inherits = FALSE)) {
      miss <- c(miss, n)
      next
    }
    a <- get(n, envir = env)
    b <- get(n, envir = ns, inherits = FALSE)
    if (!is.function(b) ||
        !identical(deparse(body(a)), deparse(body(b))) ||
        !identical(deparse(formals(a)), deparse(formals(b)))) {
      dif <- c(dif, n)
    }
  }
  cat(p, format(packageVersion(p)), dirname(find.package(p)), ": source fns",
      length(nm), " differing", length(dif), " missing", length(miss), "\n")
  if (length(dif)) cat("  differ:", head(dif, 20), "\n")
  if (length(miss)) cat("  missing:", head(miss, 20), "\n")
}
