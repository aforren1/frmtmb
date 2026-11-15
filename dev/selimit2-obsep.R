# Lane selimit2: test-se-limit.R's separated bernoulli (seed 514) under
# the BLAS this R runs: the stop, the lost reasons, the fit's warnings and
# fitted()'s warning. Under OpenBLAS 0.3.26 and 0.3.32 the fitted()
# warning took item 5423's "nocurv" words.
#   <R with OpenBLAS>/Rscript dev/selimit2-obsep.R <base|lane>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, dirname(find.package("frmtmb")), "| LAPACK", La_library(),
    "\n")
cap <- function(expr) {
  w <- character(); m <- character()
  v <- withCallingHandlers(expr,
    warning = function(z) {w <<- c(w, conditionMessage(z))
      invokeRestart("muffleWarning")},
    message = function(z) {m <<- c(m, conditionMessage(z))
      invokeRestart("muffleMessage")})
  list(v = v, w = w, m = m)
}
set.seed(514)
ds <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
ds$yb <- as.integer(ds$z > 0)
r <- cap(frm(yb ~ z + x, family = bernoulli(), data = ds))
fs <- r$v
cat("code", fs$opt$convergence, fs$opt$message, "| par",
    signif(fs$opt$par, 5), "\n")
lost <- ns$sdr_of(fs)$se_lost
cat("lost:", paste(names(lost), lost, sep = ":", collapse = ", "), "\n")
li <- ns$se_limit_info(fs)
cat("limit kinds:", paste(names(li$kind), li$kind, sep = ":",
                          collapse = ", "), "\n")
for (x in r$w) cat("FIT WARN:", gsub("\n", " ", x), "\n")
for (x in r$m) cat("FIT MSG:", gsub("\n", " ", x), "\n")
rs <- cap(stats::fitted(fs))
cat("fitted NaN", sum(is.nan(rs$v[, "Est.Error"])), "of", nrow(ds), "\n")
for (x in rs$w) cat("FITTED WARN:", gsub("\n", " ", x), "\n")
