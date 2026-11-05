# Diagnosis 1: the D3 repro on the base build (rellib-r16), verbose, so
# each round's move, re-anchored gradient and ESS are visible.
# Seeds 124 (binomial) and 69 (poisson), as dev/equiv-repro-importance.R.
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
try_fit <- function(lab, expr) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  warning:", substr(conditionMessage(w), 1, 300), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("  ", conditionMessage(m))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(r, "error")) {
    cat(lab, ": ERROR ", conditionMessage(r), "\n", sep = "")
  } else {
    cat(lab, ": logLik ", format(as.numeric(logLik(r)), digits = 10),
        if (!is.null(r$importance)) paste0(" mcse ", format(r$importance$mcse,
                                                            digits = 3)),
        " par ", paste(format(r$opt$par, digits = 6), collapse = " "),
        "\n", sep = "")
  }
  cat(sprintf("  elapsed %.2f s\n", el))
  invisible(r)
}
db <- E$eq_dat("binomial", seed = 124)
fb <- function(...) frm(bf(y | trials(ntr) ~ x + (1 | g)), family = binomial(),
                        data = db, ...)
try_fit("binomial Laplace", fb())
try_fit("binomial quadrature", fb(quadrature = TRUE))
for (k in c(200L, 512L, 1000L)) {
  try_fit(paste("binomial importance", k), fb(importance = k, verbose = TRUE))
}
dp <- E$eq_dat("poisson", seed = 69)
fp <- function(...) frm(bf(y ~ x + (1 + x | g)), family = poisson(),
                        data = dp, ...)
try_fit("poisson Laplace", fp())
for (k in c(200L, 512L)) {
  try_fit(paste("poisson importance", k), fp(importance = k, verbose = TRUE))
}
