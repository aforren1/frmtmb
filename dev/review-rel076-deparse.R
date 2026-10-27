# Review: rellib-r14 against the release tree. For each package, every
# top-level `name <- function(...)` of R/ is deparsed and compared with
# the installed namespace's object of that name; also lists top-level
# names defined twice in a package's R/, NAMESPACE exports missing from
# the installed namespace, and installed exports not in the tree's
# NAMESPACE.
#   Rscript dev/review-rel076-deparse.R > dev/review-rel076-out/deparse.txt
lib <- "C:/Users/adf44/source/r/rellib-r14"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
rel <- "C:/Users/adf44/source/r/frmtmb-wt-rel076"
pkgs <- c(frmtmb = rel,
          setNames(file.path(rel, "extensions", c("frmtmb.survival",
            "frmtmb.spatial", "frmtmb.sample", "frmtmb.latent", "frmtmb.eam",
            "frmtmb.coupling", "frmtmb.ode", "frmtmb.learn")),
            c("frmtmb.survival", "frmtmb.spatial", "frmtmb.sample",
              "frmtmb.latent", "frmtmb.eam", "frmtmb.coupling", "frmtmb.ode",
              "frmtmb.learn")))
norm <- function(f) {
  f <- utils::removeSource(f)
  environment(f) <- globalenv()
  paste(deparse(f, control = c("keepNA", "keepInteger", "niceNames",
                               "showAttributes")), collapse = "\n")
}
for (p in names(pkgs)) {
  d <- pkgs[[p]]
  ns <- asNamespace(p)
  desc <- read.dcf(file.path(d, "DESCRIPTION"))
  inst <- utils::packageVersion(p, lib.loc = lib)
  cat("==", p, "tree", desc[, "Version"], "installed", format(inst), "\n")
  files <- list.files(file.path(d, "R"), "[.][Rr]$", full.names = TRUE)
  defs <- character(0); nf <- 0; ndiff <- 0; nmiss <- 0
  for (f in files) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && (identical(e[[1]], as.name("<-")) ||
                         identical(e[[1]], as.name("="))) &&
          is.name(e[[2]])) {
        nm <- as.character(e[[2]])
        defs <- c(defs, nm)
        rhs <- e[[3]]
        if (is.call(rhs) && identical(rhs[[1]], as.name("function"))) {
          nf <- nf + 1
          src <- eval(rhs, globalenv())
          if (!exists(nm, envir = ns, inherits = FALSE)) {
            nmiss <- nmiss + 1
            cat("  MISSING in installed:", nm, "(", basename(f), ")\n")
            next
          }
          obj <- get(nm, envir = ns, inherits = FALSE)
          if (!is.function(obj) || !identical(norm(src), norm(obj))) {
            ndiff <- ndiff + 1
            cat("  DIFFERS:", nm, "(", basename(f), ")\n")
          }
        }
      }
    }
  }
  dup <- unique(defs[duplicated(defs)])
  cat("  functions compared:", nf, " differ:", ndiff, " missing:", nmiss,
      " names defined twice:", length(dup),
      if (length(dup)) paste0("(", paste(dup, collapse = ", "), ")") else "",
      "\n")
  nsf <- parseNamespaceFile(basename(d), dirname(d))
  ex_tree <- sort(nsf$exports)
  ex_inst <- sort(getNamespaceExports(ns))
  cat("  exports tree:", length(ex_tree), " installed:", length(ex_inst),
      " only tree:", paste(setdiff(ex_tree, ex_inst), collapse = ","),
      " only installed:", paste(setdiff(ex_inst, ex_tree), collapse = ","),
      "\n")
  s3t <- nrow(nsf$S3methods)
  s3i <- nrow(getNamespaceInfo(ns, "S3methods"))
  cat("  S3 methods tree:", s3t, " installed:", s3i, "\n")
}
