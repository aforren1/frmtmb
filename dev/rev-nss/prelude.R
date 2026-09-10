# Reviewer prelude for lane `nss`. LANE is the reviewer's own private
# library holding the LANE build of frmtmb.ode; REF is the round's
# shared reference build of the base commit, READ-ONLY.
LANE <- "C:/Users/adf44/source/r/rev-nss-lib"
REF  <- "C:/Users/adf44/source/r/rellib-0552"
PIN  <- "C:/Users/adf44/source/r/pinlib"
USR  <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
NONE <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
               "c--Users-adf44-source-r-frmtmb/",
               "529b6e73-d28f-46aa-a279-7dbeeb58fd4f/",
               "scratchpad/nonelib")
NOUNDONE <- sub("nonelib$", "nulib", NONE)
arm <- Sys.getenv("REV_ARM", "lane")
.libPaths(switch(arm,
                 ref  = c(REF, PIN, USR),
                 none = c(NONE, REF, PIN, USR),
                 noundone = c(NOUNDONE, REF, PIN, USR),
                 c(LANE, REF, PIN, USR)))
rev_env <- function() {
  cat("arm:", arm, "\n")
  for (p in c("frmtmb", "frmtmb.ode", "RTMB", "RTMBode", "deSolve",
              "StanHeaders")) {
    v <- tryCatch(format(packageVersion(p)), error = function(e) "-")
    w <- tryCatch(dirname(find.package(p)), error = function(e) "-")
    cat(sprintf("  %-12s %-8s %s\n", p, v, w))
  }
}
