# Reviewer, punch round 2: the library stack for the final recheck.
#
# A THIRD private library, built fresh from the worktree after the
# lane's second punch, so nothing here reads the lane's library or
# first round's `rev-latent-lib`.

REV_LIB <- "C:/Users/adf44/source/r/rev-latent-lib3"
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
