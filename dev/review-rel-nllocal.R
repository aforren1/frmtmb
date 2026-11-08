# Release review 0.80.0: can nl_local_message() (which calls logLik())
# be reached under REML with a prior on the fixed effects, where logLik()
# stops? A body that is constant at its start and a local maximum there
# (`a + b^2 * x * (1 - b^2)` on data with no x effect), without a prior,
# with a prior, and each under REML = TRUE. nl_local_message() is traced:
# a line says how many bodies it was handed.
#
#   Rscript dev/review-rel-nllocal.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
trace("nl_local_message", quote(cat("  TRACE nl_local_message bodies:",
                                    length(bodies), "\n")),
      where = ns, print = FALSE)
run <- function(lbl, expr) {
  cs <- character()
  r <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      cs <<- c(cs, paste0("W[", class(w)[1L], "] ",
                          substr(conditionMessage(w), 1, 90)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cs <<- c(cs, paste0("M[", class(m)[1L], "] ",
                          substr(conditionMessage(m), 1, 90)))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      cs <<- c(cs, paste0("ERROR[", class(e)[1L], "] ",
                          substr(conditionMessage(e), 1, 160)))
      NULL
    })
  cat(lbl, if (is.null(r)) "no fit" else
    sprintf("penalized %.4f", as.numeric(logLik_penalized(r))), "\n")
  if (length(cs)) cat(paste0("  ", cs), sep = "\n")
}
set.seed(4)
n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 1 + rnorm(n)
fm <- bf(y ~ a + b^2 * x * (1 - b^2), a ~ 1, b ~ 1, nl = TRUE)
pr <- prior(normal(0, 10), nlpar = "a")
run("ML, no prior   ", frm(fm, data = d))
run("ML, prior on a ", frm(fm, data = d, prior = pr))
run("REML, no prior ", frm(fm, data = d, REML = TRUE))
run("REML, prior on a", frm(fm, data = d, prior = pr, REML = TRUE))
cat("NLLOCAL DONE\n")
