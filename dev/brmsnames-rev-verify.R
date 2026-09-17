## Reviewer: does the lane's installed build match the worktree source,
## function by function? Rscript dev/brmsnames-rev-verify.R
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
chk <- function(pkg, dir) {
  ns <- asNamespace(pkg)
  env <- new.env(parent = ns)
  files <- sort(list.files(file.path(dir, "R"), "[.]R$", full.names = TRUE))
  for (f in files) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      # only top-level assignments of functions; skip side effects
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
          is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
        eval(e, env)
    }
  }
  nms <- ls(env, all.names = TRUE)
  bad <- character(0); miss <- character(0)
  for (n in nms) {
    if (!exists(n, envir = ns, inherits = FALSE)) { miss <- c(miss, n); next }
    a <- get(n, envir = ns)
    if (bindingIsActive(n, ns)) next
    b <- get(n, envir = env)
    if (!is.function(a)) { bad <- c(bad, n); next }
    da <- deparse(removeSource(a)); db <- deparse(removeSource(b))
    if (!identical(da, db)) bad <- c(bad, n)
  }
  cat(pkg, ": source functions", length(nms), " differ", length(bad),
      " missing", length(miss), "\n")
  if (length(bad)) print(bad)
  if (length(miss)) print(miss)
}
chk("frmtmb", ".")
chk("frmtmb.sample", "extensions/frmtmb.sample")
chk("frmtmb.learn", "extensions/frmtmb.learn")
if (dir.exists("C:/Users/adf44/source/r/brmsnames-lib/frmtmb.latent")) chk("frmtmb.latent", "extensions/frmtmb.latent")
cat("versions:", format(packageVersion("frmtmb")),
    format(packageVersion("frmtmb.sample")), "\n")
cat("lib frmtmb:", find.package("frmtmb"), "\n")
