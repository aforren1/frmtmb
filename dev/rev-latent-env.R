# Reviewer `latent`: the library stack every rev-latent-* script sources.
#
# The reviewer builds its OWN copy of frmtmb.latent from the worktree,
# because the lane's own library was at one point stale (its installer
# gated on PowerShell's `$?`), and a review that reads the lane's
# library is reviewing whatever happened to be installed there.
# Everything else comes from the round's shared reference build.

REV_LIB <- "C:/Users/adf44/source/r/rev-latent-lib"
dir.create(REV_LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(REV_LIB,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

rev_env_report <- function() {
  cat("libPaths:\n")
  cat(paste0("  ", .libPaths(), collapse = "\n"), "\n")
  for (p in c("frmtmb", "frmtmb.latent", "RTMB", "TMB",
              "depmixS4", "hmmTMB", "poLCA", "StanHeaders")) {
    loc <- find.package(p, quiet = TRUE)
    cat(sprintf("  %-14s %-10s %s\n", p,
                if (length(loc)) as.character(utils::packageVersion(p))
                else "MISSING",
                if (length(loc)) loc[1L] else ""))
  }
  invisible(TRUE)
}
