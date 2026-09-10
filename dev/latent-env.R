# Lane `latent`: the library stack every script in this lane sources.
#
# Order matters and is the one dev/lane-rules.md pins: this lane's own
# library first (it holds the change), then the round's shared reference
# build of the base commit, then the StanHeaders pin, then the user
# library as a read-only fallback. Nothing here ever installs.

LANE_LIB <- "C:/Users/adf44/source/r/latent-lib"
dir.create(LANE_LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LANE_LIB,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

lane_env_report <- function() {
  cat("libPaths:\n"); cat(paste0("  ", .libPaths(), collapse = "\n"), "\n")
  for (p in c("frmtmb", "frmtmb.latent", "RTMB", "TMB",
              "depmixS4", "hmmTMB", "poLCA")) {
    loc <- find.package(p, quiet = TRUE)
    cat(sprintf("  %-14s %-10s %s\n", p,
                if (length(loc)) as.character(utils::packageVersion(p))
                else "MISSING",
                if (length(loc)) loc[1L] else ""))
  }
  invisible(TRUE)
}
