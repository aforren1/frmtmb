# Roxygenise and install core, then every extension, from this worktree
# into the lane library. Core first: an extension's roxygenise loads the
# core it will run against.
#
# Usage: Rscript dev/superassign-install.R [pkg ...]   (default: all)
lib <- "C:/Users/adf44/source/r/wt-superassign-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(lib, base, user))
root <- "C:/Users/adf44/source/r/frmtmb-wt-superassign"
all <- c("frmtmb", paste0("frmtmb.", c("eam", "sample", "coupling",
  "learn", "latent", "ode", "spatial", "survival")))
want <- commandArgs(trailingOnly = TRUE)
if (!length(want)) want <- all
path_of <- function(p) {
  if (p == "frmtmb") root else file.path(root, "extensions", p)
}
lp <- sprintf(".libPaths(c('%s','%s','%s'))", lib, base, user)
for (p in intersect(all, want)) {
  dir <- path_of(p)
  cat("== roxygenise", p, "\n")
  out <- system2(file.path(R.home("bin"), "Rscript"),
                 c("-e", shQuote(sprintf("%s; roxygen2::roxygenise('%s')",
                                         lp, dir))),
                 stdout = TRUE, stderr = TRUE)
  cat(tail(out, 3), sep = "\n")
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) stop("roxygenise failed: ", p)
  cat("== install", p, "\n")
  out <- system2(file.path(R.home("bin"), "R"),
                 c("CMD", "INSTALL", paste0("--library=", lib),
                   "--no-multiarch", shQuote(dir)),
                 stdout = TRUE, stderr = TRUE)
  cat(tail(out, 2), sep = "\n")
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) stop("install failed: ", p)
}
for (p in intersect(all, want)) {
  cat(p, as.character(packageVersion(p, lib.loc = lib)), "\n")
}
