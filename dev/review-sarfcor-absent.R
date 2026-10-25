# Reviewer: a saved sar() fit read where frmtmb.spatial is not installed
# (absent) and where it is installed but cannot load (unloadable), with
# the helper the extensions use. Every method should refuse by the
# package's name (class frmtmb_extension_missing).
# Rscript dev/review-sarfcor-absent.R
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
skip_env <- new.env()
skip_env$testthat <- NULL
sys.source(file.path(wt, "extensions/frmtmb.spatial/tests/testthat",
                     "helper-absent-package.R"), envir = environment())
set.seed(3)
g <- expand.grid(r = 1:5, c = 1:5)
n <- nrow(g)
A <- outer(1:n, 1:n, function(i, j)
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
W <- A / rowSums(A)
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.5 * W, 1 + 0.5 * d$x + rnorm(n)))
fit <- frm(y ~ x + sar(W), data = d, data2 = list(W = W))
rds <- tempfile(fileext = ".rds")
saveRDS(fit, rds)
code <- c(
  "suppressMessages(library(frmtmb))",
  "cat('spatial installed:', nzchar(system.file(package = 'frmtmb.spatial')), '\\n')",
  sprintf("fit <- readRDS('%s')", gsub("\\\\", "/", rds)),
  "for (m in c('print', 'summary', 'logLik', 'fitted', 'predict', 'residuals', 'simulate', 'autocor_matrix', 'confint_varcorr', 'variables', 'vcov', 'fixef', 'update')) {",
  "  r <- tryCatch({x <- if (m == 'update') update(fit, data = fit$data) else",
  "    if (m == 'print') utils::capture.output(print(fit)) else",
  "    if (m == 'summary') utils::capture.output(print(summary(fit))) else",
  "    get(m)(fit); 'ANSWERED'},",
  "    error = function(e) paste(if (inherits(e, 'frmtmb_extension_missing')) 'MISSING' else 'OTHER',",
  "      substr(conditionMessage(e), 1, 140)))",
  "  cat(sprintf('%-16s %s\\n', m, r))",
  "}")
for (arm in c("unloadable", "absent")) {
  libs <- if (arm == "absent") {
    lib <- tempfile("absent-lib"); dir.create(lib)
    file.copy(find.package("frmtmb"), lib, recursive = TRUE)
    c(lib, Filter(function(l) !dir.exists(file.path(l, "frmtmb.spatial")),
                  .libPaths()))
  } else unloadable_package_libs("frmtmb.spatial")
  cat("==== arm:", arm, "\n")
  cat(child_with_libs(libs, code), sep = "\n")
}
