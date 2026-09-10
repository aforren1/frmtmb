# Shared by every script in this directory. Sourced, not run.
#
# The library order is the one dev/lane-rules.md fixes: this lane's own
# private library first, then the pinned StanHeaders, then the round's
# shared reference build of the base commit, then the user library as a
# last resort. RLDDM_LIB selects the arm, so a script can be run against
# the reference build by setting it.

lane_lib <- Sys.getenv("RLDDM_LIB", "C:/Users/adf44/source/r/rlddm-lib")
.libPaths(c(lane_lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

cat("lib      :", lane_lib, "\n")
cat("frmtmb   :", format(packageVersion("frmtmb")), "at",
    dirname(system.file(package = "frmtmb")), "\n")
cat("eam      :", format(packageVersion("frmtmb.eam")), "at",
    dirname(system.file(package = "frmtmb.eam")), "\n")
cat("learn    :", format(packageVersion("frmtmb.learn")), "at",
    dirname(system.file(package = "frmtmb.learn")), "\n")
