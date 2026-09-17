# What the shared libraries and CRAN hold today, before anything is
# installed. Read-only: installs nothing.
pk <- c("tmbstan", "rstan", "StanHeaders", "RcppParallel", "BH", "TMB",
        "RTMB", "RcppEigen", "Rcpp", "rstantools", "inline", "QuickJSR",
        "brms", "loo", "frmtmb", "frmtmb.sample")
libs <- c(pin = "C:/Users/adf44/source/r/pinlib",
          user = "C:/Users/adf44/AppData/Local/R/win-library/4.6")
for (nm in names(libs)) {
  ip <- installed.packages(lib.loc = libs[[nm]], fields = "Built")
  p <- intersect(pk, rownames(ip))
  cat("==", nm, libs[[nm]], "\n")
  print(ip[p, c("Version", "Built"), drop = FALSE])
}
repo <- "https://cloud.r-project.org"
src <- available.packages(repos = repo, type = "source")
bin <- available.packages(repos = repo, type = "win.binary")
q <- c("tmbstan", "StanHeaders", "rstan", "TMB", "RTMB", "RcppParallel",
       "BH", "RcppEigen")
print(cbind(source = src[q, "Version"], win.binary = bin[q, "Version"]))
