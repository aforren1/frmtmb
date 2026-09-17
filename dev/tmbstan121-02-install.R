# Build one arm's private library. Usage:
#   Rscript tmbstan121-02-install.R <arm>
# arm is one of common, A, Bsrc, Bbin, C120, C120A.
#
# Every arm lives under ROOT, which is this lane's private library. The
# StanHeaders a tmbstan build sees is decided by .libPaths() at install
# time, because tmbstan's configure runs rstan::stanc() to regenerate
# model.hpp, and stanc comes from whichever StanHeaders loads first.
ROOT <- "C:/Users/adf44/source/r/tmbstan121-lib"
PIN <- "C:/Users/adf44/source/r/pinlib"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
SRC <- file.path(ROOT, "src")
WT <- "C:/Users/adf44/source/r/frmtmb-wt-tmbstan121"

arm <- commandArgs(trailingOnly = TRUE)[1]
LIB <- file.path(ROOT, arm)
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)

paths <- switch(arm,
  common = c(LIB, PIN, USER),
  A      = c(LIB, PIN, USER),
  C120A  = c(LIB, PIN, USER),
  Bsrc   = c(LIB, USER),
  Bbin   = c(LIB, USER),
  C120   = c(LIB, USER),
  stop("unknown arm ", arm))
.libPaths(paths)
Sys.setenv(R_LIBS = paste(paths, collapse = ";"))
cat("arm", arm, "libPaths:\n")
print(.libPaths())
cat("StanHeaders at install:", format(packageVersion("StanHeaders")),
    "from", find.package("StanHeaders"), "\n")

inst <- function(pkg, type = "source") {
  install.packages(pkg, lib = LIB, repos = NULL, type = type,
                   INSTALL_opts = "--no-multiarch")
}

if (arm == "common") {
  for (p in c(".", "extensions/frmtmb.eam", "extensions/frmtmb.latent",
              "extensions/frmtmb.ode", "extensions/frmtmb.sample")) {
    inst(normalizePath(file.path(WT, p)))
  }
} else if (arm %in% c("A", "Bsrc")) {
  inst(file.path(SRC, "tmbstan_1.2.1.tar.gz"))
} else if (arm == "Bbin") {
  inst(file.path(SRC, "tmbstan_1.2.1.zip"), type = "win.binary")
} else {
  inst(file.path(SRC, "tmbstan_1.2.0.tar.gz"))
}

ip <- installed.packages(lib.loc = LIB, fields = "Built")
print(ip[, c("Version", "Built"), drop = FALSE])
hpp <- file.path(LIB, "tmbstan", "model.hpp")
if (file.exists(hpp)) {
  x <- readLines(hpp, warn = FALSE)
  cat("model.hpp lines:", length(x),
      " log_prob_impl overloads:", sum(grepl("^  log_prob_impl[(]", x)),
      " patched custom_func(y):",
      sum(grepl("custom_func::custom_func(y)", x, fixed = TRUE)),
      " UNPATCHED std_normal_lpdf<propto__>(y):",
      sum(grepl("std_normal_lpdf<propto__>(y)", x, fixed = TRUE)), "\n")
}
