# Release review of 0.74.0, item 1: does rellib-r12 hold the tree?
# For each package: every top-level `name <- function(...)` in the
# tree's R/ is deparsed and compared with the installed namespace's
# object; names on one side only are listed; the installed NAMESPACE,
# DESCRIPTION version and help topics are compared with the tree's.
#   Rscript dev/review-rel074-libmatch.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
pk <- c(frmtmb = ".", setNames(file.path("extensions",
  c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
    "frmtmb.ode", "frmtmb.sample", "frmtmb.survival", "frmtmb.spatial")),
  c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
    "frmtmb.ode", "frmtmb.sample", "frmtmb.survival", "frmtmb.spatial")))
strip <- function(f) {
  if (is.function(f)) {
    attributes(f) <- NULL
    f <- removeSource(f)
  }
  paste(deparse(f, control = c("keepNA", "keepInteger", "niceNames",
                                "showAttributes")), collapse = "\n")
}
for (p in names(pk)) {
  dir <- file.path(root, pk[[p]])
  ns <- asNamespace(p)
  env <- new.env(parent = ns)
  defs <- character()
  dup <- character()
  for (f in list.files(file.path(dir, "R"), "[.][Rr]$", full.names = TRUE)) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) {
        nm <- as.character(e[[2]])
        if (nm %in% defs) dup <- c(dup, paste0(nm, " (", basename(f), ")"))
        defs <- c(defs, nm)
        assign(nm, eval(e[[3]], env), envir = env)
      }
    }
  }
  nsfun <- Filter(function(n) is.function(get(n, ns)), ls(ns, all.names = TRUE))
  diff <- character()
  for (n in unique(defs)) {
    if (!exists(n, ns, inherits = FALSE)) {
      diff <- c(diff, paste0("MISSING-IN-LIB ", n))
      next
    }
    a <- strip(get(n, env)); b <- strip(get(n, ns))
    if (!identical(a, b)) diff <- c(diff, paste0("DIFF ", n))
  }
  extra <- setdiff(nsfun, defs)
  cat("==", p, as.character(packageVersion(p)), "tree DESCRIPTION",
      read.dcf(file.path(dir, "DESCRIPTION"), "Version"), "\n")
  cat("functions in tree:", length(unique(defs)), "compared; differ:",
      length(diff), "; duplicate top-level definitions:", length(dup), "\n")
  if (length(diff)) cat(paste0("  ", diff), sep = "\n")
  if (length(dup)) cat(paste0("  DUP ", dup), sep = "\n")
  cat("functions in lib only (not a plain top-level def):", length(extra),
      if (length(extra)) paste(head(extra, 30), collapse = " "), "\n")
  # NAMESPACE
  a <- readLines(file.path(dir, "NAMESPACE"))
  b <- readLines(file.path(find.package(p), "NAMESPACE"))
  cat("NAMESPACE identical:", identical(a, b), "\n")
  # DESCRIPTION fields
  d1 <- read.dcf(file.path(dir, "DESCRIPTION"))
  d2 <- read.dcf(file.path(find.package(p), "DESCRIPTION"))
  for (fl in c("Version", "Depends", "Imports", "Suggests", "LinkingTo")) {
    x <- if (fl %in% colnames(d1)) gsub("\\s+", " ", d1[, fl]) else NA
    y <- if (fl %in% colnames(d2)) gsub("\\s+", " ", d2[, fl]) else NA
    if (!identical(x, y)) cat("DESCRIPTION", fl, "differs:\n  tree:", x,
                              "\n  lib: ", y, "\n")
  }
  # help topics
  rd <- sub("[.]Rd$", "", list.files(file.path(dir, "man"), "[.]Rd$"))
  idx <- readRDS(file.path(find.package(p), "help", "aliases.rds"))
  libtop <- unique(unname(idx))
  cat("Rd files in tree:", length(rd), "; help files in lib:",
      length(libtop), "; only tree:",
      paste(setdiff(rd, libtop), collapse = " "), "; only lib:",
      paste(setdiff(libtop, rd), collapse = " "), "\n")
  # vignettes
  v1 <- list.files(file.path(dir, "vignettes"), "[.]Rmd$")
  v2 <- if (dir.exists(file.path(find.package(p), "doc")))
    list.files(file.path(find.package(p), "doc"), "[.]Rmd$") else character()
  cat("vignettes tree:", paste(v1, collapse = " "), "| lib doc:",
      paste(v2, collapse = " "), "\n")
}
