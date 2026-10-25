# The names of every fit of the core and frmtmb.sample suites on the
# merged build (dev/rel075-equality.sh suite, recorded by a copy of
# dev/brmsnames-1024-suitenames.R) against the same record on the
# brmsnames lane's build (its dev/brmsnames-1024-suite/lane/, read in
# the lane's tree). Fits pair by test file and call order. A file whose
# names differ is listed with the lanes that changed it, from the
# merge's status lists (dev/rel075-log/merge-work/<lane>.status and
# copied-<lane>.lst): a difference is accounted for only in a file
# another lane changed or added.
#   Rscript dev/rel075-eq-suitenames.R
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel075")
M <- "dev/rel075-log/eq/bn-suite/lane"
L <- "C:/Users/adf44/source/r/frmtmb-wt-brmsnames/dev/brmsnames-1024-suite/lane"
W <- "dev/rel075-log/merge-work"
touched <- list()
for (l in c("smallfix", "bhaz", "brmsnames", "sarfcor")) {
  st <- read.delim(file.path(W, paste0(l, ".status")), header = FALSE)
  touched[[l]] <- basename(c(st[[2]], readLines(file.path(
    "dev/rel075-log", paste0("copied-", l, ".lst")))))
}
fields <- c("internal", "brms", "variables", "outer", "prior")
for (pkg in c("frmtmb", "frmtmb.sample")) {
  fm <- list.files(file.path(M, pkg), "[.]rds$")
  fl <- list.files(file.path(L, pkg), "[.]rds$")
  both <- intersect(fm, fl)
  cat("==", pkg, ": files merged", length(fm), "| lane", length(fl),
      "| both", length(both), "| only merged:",
      paste(setdiff(fm, fl), collapse = ", "), "| only lane:",
      paste(setdiff(fl, fm), collapse = ", "), "\n")
  nfit <- 0L; nsame <- 0L; bad <- character(0)
  for (f in both) {
    a <- readRDS(file.path(M, pkg, f)); b <- readRDS(file.path(L, pkg, f))
    k <- min(length(a$fits), length(b$fits))
    same <- vapply(seq_len(k), function(i) {
      identical(a$fits[[i]][fields], b$fits[[i]][fields])
    }, NA)
    nfit <- nfit + k; nsame <- nsame + sum(same)
    if (!all(same) || length(a$fits) != length(b$fits) ||
          !identical(a$counts, b$counts)) {
      tf <- sub("[.]rds$", ".R", f)
      by <- names(Filter(function(v) tf %in% v, touched))
      bad <- c(bad, sprintf(
        "  %s: fits merged %d lane %d, %d of %d paired fits differ; counts merged [%s] lane [%s]; test file changed by: %s",
        f, length(a$fits), length(b$fits), sum(!same), k, a$counts,
        b$counts, if (length(by)) paste(by, collapse = ", ") else "NONE"))
    }
  }
  cat("paired fits", nfit, "| identical in all five name fields", nsame,
      "\nfiles that differ:", length(bad), "\n")
  if (length(bad)) cat(bad, sep = "\n")
}
