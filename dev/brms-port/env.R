# Where the audit reads and writes, and which frmtmb build it measures.
#
# Why environment variables: the audit is rerun against more than one
# build (a release library, and older ones to date a regression), and
# each run must land in its own directory. Up to 0.34.0 every path here
# was a constant naming one worktree, so the scripts broke when that
# worktree was removed.
#
#   PORT_LIB  libraries to put before the user library, ";"-separated.
#             Unset: measure whatever frmtmb .libPaths() already finds.
#   PORT_OUT  directory for the per-vignette RDS and logs and for the
#             merged tables. Unset: dev/brms-port-out, which
#             dev/.gitignore keeps out of the repository.
#
# The caller sets HERE (this directory) before sourcing this file, with
#   HERE <- local({
#     a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
#     normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
#   })

PORT_ROOT <- dirname(dirname(HERE))
PORT_OUT <- Sys.getenv("PORT_OUT",
                       unset = file.path(PORT_ROOT, "dev", "brms-port-out"))
dir.create(PORT_OUT, showWarnings = FALSE, recursive = TRUE)
local({
  lib <- Sys.getenv("PORT_LIB", unset = "")
  if (nzchar(lib)) {
    .libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
                "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  }
})

# The build under test, recorded in every log so a number cannot be
# read against the wrong frmtmb.
port_build <- function() {
  sprintf("frmtmb %s from %s", utils::packageVersion("frmtmb"),
          dirname(system.file(package = "frmtmb")))
}
