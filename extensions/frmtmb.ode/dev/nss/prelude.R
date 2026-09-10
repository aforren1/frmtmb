# Shared prelude for lane `nss` scripts.
#
# LIB is the lane's private library and is the ONLY place anything is
# installed. REF is the round's shared reference build of the base
# commit and is READ-ONLY; a script that wants the BEFORE behaviour
# sets NSS_ARM=ref before sourcing this.

LIB <- "C:/Users/adf44/source/r/nsslib"
REF <- "C:/Users/adf44/source/r/rellib-0552"
PIN <- "C:/Users/adf44/source/r/pinlib"
USR <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"

arm <- Sys.getenv("NSS_ARM", "lane")
.libPaths(if (identical(arm, "ref")) c(REF, PIN, USR) else
          c(LIB, REF, PIN, USR))

nss_report_env <- function() {
  cat("arm:", arm, "\n")
  for (p in c("frmtmb", "frmtmb.ode", "RTMB", "RTMBode", "deSolve")) {
    v <- tryCatch(format(packageVersion(p)), error = function(e) "-")
    w <- tryCatch(dirname(find.package(p)), error = function(e) "-")
    cat(sprintf("  %-12s %-8s %s\n", p, v, w))
  }
  cat("  R", R.version.string, "\n")
}
