# Reviewer of the 0.70.0 consolidation: the ported row brmsfit-methods:326
# (all(fitted(fit1, dpar = "sigma") > 0)) on fixture 1, per build: the
# optimizer's code and message, the objective, the smallest log sd, and
# whether the row holds.
#   Rscript dev/relrev070-row326.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
T <- "C:/Users/adf44/source/r/frmtmb-wt-release"
cat("frmtmb", as.character(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
# fixture 1 as the ported suite builds it (helper-brms-suite.R)
env <- new.env()
sys.source(file.path(T, "tests/testthat/helper-brms-suite.R"), envir = env)
fx <- get("brms_fixture", envir = env)
fit1 <- suppressWarnings(fx(1))
fi <- suppressWarnings(fitted(fit1, dpar = "sigma"))
th <- fit1$opt$par[names(fit1$opt$par) == "theta"]
cat("code", fit1$opt$convergence, "|", fit1$opt$message, "| objective",
    sprintf("%.10f", fit1$opt$objective), "| min theta",
    format(min(th), digits = 5), "\n")
cat("fitted sigma columns:", paste(colnames(fi), collapse = ","),
    "| NaN per column:", paste(colSums(is.nan(as.matrix(fi))), collapse = ","),
    "| all(fi > 0):", all(fi > 0), "| Estimate > 0:",
    all(fi[, 1] > 0), "\n")
