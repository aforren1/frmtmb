# Reviewer: does the lane library hold the worktree's source?
# Every function sourced from R/ (collate order) is compared by deparse
# with the installed namespace's binding; also NAMESPACE exports.
LIB <- "C:/Users/adf44/source/r/wt-formrobust-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust"
chk <- function(pkg, dir) {
  ns <- asNamespace(pkg)
  cat(pkg, "from", find.package(pkg), "\n")
  desc <- read.dcf(file.path(dir, "DESCRIPTION"))
  files <- if ("Collate" %in% colnames(desc)) {
    scan(text = gsub("\n", " ", desc[, "Collate"]), what = "",
         quiet = TRUE)
  } else sort(list.files(file.path(dir, "R"), pattern = "[.][Rr]$"))
  e <- new.env(parent = ns)
  for (f in files) sys.source(file.path(dir, "R", f), envir = e,
                               keep.source = FALSE)
  nms <- ls(e, all.names = TRUE)
  bad <- character(0); nf <- 0L; missing <- character(0)
  for (nm in nms) {
    x <- get(nm, envir = e)
    if (!is.function(x)) next
    nf <- nf + 1L
    if (!exists(nm, envir = ns, inherits = FALSE)) {
      missing <- c(missing, nm); next
    }
    if (bindingIsActive(nm, ns)) next
    y <- get(nm, envir = ns, inherits = FALSE)
    if (!identical(deparse(x), deparse(y))) bad <- c(bad, nm)
  }
  cat(pkg, ": functions sourced", nf, " differing", length(bad),
      " missing", length(missing), "\n")
  if (length(bad)) cat("  DIFF:", bad, "\n")
  if (length(missing)) cat("  MISSING:", missing, "\n")
  inst_ns <- readLines(file.path(find.package(pkg), "NAMESPACE"))
  src_ns <- readLines(file.path(dir, "NAMESPACE"))
  cat(pkg, ": NAMESPACE identical:", identical(inst_ns, src_ns), "\n")
}
chk("frmtmb", WT)
chk("frmtmb.sample", file.path(WT, "extensions/frmtmb.sample"))
