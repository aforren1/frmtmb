# Lane `learnhier` (item 2.2): the library stack every script here
# sources.
#
# Order is the one dev/lane-rules.md pins: this lane's own library first
# (it holds the change), then the round's shared reference build of the
# base commit, then the StanHeaders pin, then the user library as a
# read-only fallback. Nothing here ever installs.

LANE_LIB <- "C:/Users/adf44/source/r/learnhier-lib"
dir.create(LANE_LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LANE_LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

# The Stan cache in this worktree. Without it every program in the
# identity tier recompiles, and rstan 2.32.7 cannot compile against the
# user library's StanHeaders 2.39.1, so a fresh compile is also the one
# place the pin can fail quietly.
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE =
               "C:/Users/adf44/source/r/frmtmb-wt-learnhier/dev/stan-cache")
}

lane_env_report <- function() {
  cat("libPaths:\n")
  cat(paste0("  ", .libPaths(), collapse = "\n"), "\n")
  for (p in c("frmtmb", "frmtmb.learn", "frmtmb.eam", "RTMB", "TMB",
              "rstan", "StanHeaders", "RWiener")) {
    loc <- find.package(p, quiet = TRUE)
    cat(sprintf("  %-14s %-10s %s\n", p,
                if (length(loc)) as.character(utils::packageVersion(p))
                else "MISSING",
                if (length(loc)) loc[1L] else ""))
  }
  cat("  FRMTMB_STAN_CACHE ", Sys.getenv("FRMTMB_STAN_CACHE"), "\n")
  invisible(TRUE)
}
