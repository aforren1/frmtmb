# Reviewer of lane setier: copy of dev/nanse-rev2-ridge-re.R with lane = wt-setier-lib over
# rellib-r6, base = rellib-r6 (merge arm = rellib-r6 alone, i.e. base).
# Reviewer, punch round 1: an exact ridge c + dd on a model WITH random
# effects (finite-difference Hessian). Design B of nanse-rev2-b1.R at
# cor 0.999 reported no lost SE. What SEs do c and dd get, on the lane
# and on base, and over several seeds?
#   Rscript dev/nanse-rev2-ridge-re.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
for (s in 1:12) {
  set.seed(s)
  n <- 300
  x1 <- rnorm(n)
  g <- factor(rep(1:15, 20))
  d <- data.frame(x1, g, y = 1 + 0.5 * x1 + rnorm(15, 0, 0.5)[g] + rnorm(n))
  w <- character()
  f <- withCallingHandlers(
    frm(bf(y ~ a + c + dd, a ~ 0 + x1, c ~ 1 + (1 | g), dd ~ 1, nl = TRUE),
        data = d),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  V <- suppressWarnings(vcov(f, full = TRUE))
  se <- suppressWarnings(sqrt(diag(V)))
  names(se) <- ns$outer_par_names(f)
  cat(sprintf("seed %2d code %d: SE c %.3g dd %.3g a %.3g | lost %s | warnings %d%s\n",
              s, f$opt$convergence, se[["c_(Intercept)"]],
              se[["dd_(Intercept)"]], se[["a_x1"]],
              paste(names(ns$sdr_of(f)$se_lost), collapse = ","), length(w),
              if (length(w)) paste0(" [", substr(w[1], 1, 50), "]") else ""))
}
