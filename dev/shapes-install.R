# Install this worktree's packages into the lane's PRIVATE library.
#
# The lane library holds the lane's change; the shared read-only build
# `rellib-r3` holds the base commit it is measured against. Never
# install into rellib-r3 or into the user library.

LIB <- "C:/Users/adf44/source/r/shapes-lib"
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

wt <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
args <- commandArgs(trailingOnly = TRUE)
which <- if (length(args)) args else c("core", "sample")

run <- function(path) {
  cat("== installing", path, "\n")
  st <- system2("R", c("CMD", "INSTALL", paste0("--library=", LIB),
                       shQuote(path)), stdout = TRUE, stderr = TRUE)
  cat(tail(st, 6), sep = "\n")
  invisible(NULL)
}

if ("core" %in% which) run(wt)
if ("sample" %in% which) run(file.path(wt, "extensions/frmtmb.sample"))
for (ex in setdiff(which, c("core", "sample"))) {
  run(file.path(wt, "extensions", ex))
}
cat("frmtmb", as.character(utils::packageVersion("frmtmb", lib.loc = LIB)),
    "\n")
