# Reviewer: does the worker's installed build match the worktree source?
# Compares the deparsed body of every function defined in the changed R
# files against the same name in the installed namespace.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("frmtmb lib:", dirname(system.file(package = "frmtmb")), "\n")
cat("frmtmb version:", as.character(packageVersion("frmtmb")), "\n")
cat("spline lib:", dirname(system.file(package = "frmtmb.spline")), "\n")
cat("sample lib:", dirname(system.file(package = "frmtmb.sample")), "\n")

root <- "C:/Users/adf44/source/r/frmtmb-wt-resmooth"
files <- c("R/autocor.R", "R/bootstrap.R", "R/conditional-effects.R",
           "R/interop.R", "R/predict-brms.R", "R/predict.R",
           "R/re-formula.R", "R/simulate-newdata.R")
ns <- asNamespace("frmtmb")
bad <- character(0)
nchk <- 0L
for (f in files) {
  e <- new.env()
  sys.source(file.path(root, f), envir = e, keep.source = FALSE)
  for (nm in ls(e, all.names = TRUE)) {
    src <- get(nm, envir = e)
    if (!is.function(src)) next
    if (!exists(nm, envir = ns, inherits = FALSE)) {
      bad <- c(bad, paste0(f, ":", nm, " MISSING in namespace")); next
    }
    inst <- get(nm, envir = ns, inherits = FALSE)
    nchk <- nchk + 1L
    a <- paste(deparse(body(src)), collapse = "\n")
    b <- paste(deparse(body(inst)), collapse = "\n")
    fa <- paste(deparse(formals(src)), collapse = "\n")
    fb <- paste(deparse(formals(inst)), collapse = "\n")
    if (!identical(a, b) || !identical(fa, fb)) {
      bad <- c(bad, paste0(f, ":", nm, " DIFFERS"))
    }
  }
}
cat("functions compared:", nchk, "\n")
cat("mismatches:", length(bad), "\n")
if (length(bad)) cat(paste(bad, collapse = "\n"), "\n")

# spline
e <- new.env()
sys.source(file.path(root, "extensions/frmtmb.spline/R/curve.R"),
           envir = e, keep.source = FALSE)
nss <- asNamespace("frmtmb.spline")
badd <- character(0); n2 <- 0L
for (nm in ls(e, all.names = TRUE)) {
  src <- get(nm, envir = e)
  if (!is.function(src)) next
  if (!exists(nm, envir = nss, inherits = FALSE)) {
    badd <- c(badd, paste0("curve.R:", nm, " MISSING")); next
  }
  n2 <- n2 + 1L
  if (!identical(paste(deparse(body(src)), collapse = "\n"),
                 paste(deparse(body(get(nm, envir = nss))), collapse = "\n"))) {
    badd <- c(badd, paste0("curve.R:", nm, " DIFFERS"))
  }
}
cat("spline functions compared:", n2, "mismatches:", length(badd), "\n")
if (length(badd)) cat(paste(badd, collapse = "\n"), "\n")

cat("StanHeaders:", as.character(packageVersion("StanHeaders")), "\n")
cat("rstan:", as.character(packageVersion("rstan")), "\n")
cat("brms:", as.character(packageVersion("brms")), "\n")
cat("mgcv:", as.character(packageVersion("mgcv")), "\n")
cat("base frmtmb version:",
    as.character(packageVersion("frmtmb",
      lib.loc = "C:/Users/adf44/source/r/rellib-r3")), "\n")
