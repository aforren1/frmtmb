# test-se-limit.R:200 under OpenBLAS 0.3.26 with LAPACK: the separated
# bernoulli of "a fitted value that barely moves along a limit keeps its
# SE" (seed 514) gave its NaN fitted SEs without the warning. Print the SE
# check's verdict, the null basis and its "limit" marks, the NaN count and
# every condition of fitted(), under whatever BLAS this R runs.
#
#   Rscript dev/rel082-ob-sep.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("blas:", extSoftVersion()[["BLAS"]], La_library(), "\n")
conds <- function(expr) {
  got <- character()
  v <- withCallingHandlers(expr,
    warning = function(w) {
      got <<- c(got, paste("W:", substr(conditionMessage(w), 1, 100)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      got <<- c(got, paste("M:", substr(conditionMessage(m), 1, 100)))
      invokeRestart("muffleMessage")
    })
  list(value = v, conds = got)
}
set.seed(514)
ds <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
ds$yb <- as.integer(ds$z > 0)
fs <- conds(frm(yb ~ z + x, family = bernoulli(), data = ds))
cat("fit conditions:\n"); cat(paste0("  ", fs$conds, "\n"), sep = "")
f <- fs$value
sdr <- f$sdr
cat("se_lost:", paste(names(sdr$se_lost), sdr$se_lost, sep = "=", collapse = " "),
    "\n")
nb <- sdr$se_null
cat("se_null columns:", if (is.null(nb)) 0 else NCOL(nb), " limit attr:",
    paste(attr(nb, "limit"), collapse = " "), "\n")
rs <- conds(stats::fitted(f))
cat("fitted NaN:", sum(is.nan(rs$value[, "Est.Error"])), "of",
    nrow(rs$value), "\n")
cat("fitted conditions:\n"); cat(paste0("  ", rs$conds, "\n"), sep = "")
cat("estimates:", format(fixef(f)[, "Estimate"], digits = 6), "\n")
