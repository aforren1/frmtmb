# Reviewer of the 0.69.0 consolidation: is the prediction warning of
# ported row data-helpers:7 (fixture 1, OpenBLAS) a correct
# estimability warning or a false alarm? Prints the lost directions the
# prediction test reads, their loadings, the cosine of each prediction
# row with each, the design relation of Age and s(Age)'s fixed column,
# and whether the likelihood and the linear predictor move along the
# direction (finite steps, inner problem re-solved by obj$fn).
#   Rscript dev/relrev069-dh7diag.R <lib dir>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "| R", R.home(), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat"
for (h in list.files(root, "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
ns <- asNamespace("frmtmb")
withr::local_seed(3L)
fit <- suppressWarnings(suppressMessages(env$brms_fixture(1)))
cat("code", fit$opt$convergence, "objective",
    sprintf("%.10f", fit$opt$objective), "\n")
on <- ns$outer_par_names(fit)
sdr <- suppressWarnings(ns$sdr_of(fit))
cat("se_lost:\n"); print(sdr$se_lost)
jc <- ns$get_joint_cov(fit)
cat("joint-cov null basis columns:", NCOL(jc$null), "\n")
if (NCOL(jc$null)) {
  N <- as.matrix(jc$null)
  rownames(N) <- jc$names
  for (j in seq_len(ncol(N))) {
    nz <- which(abs(N[, j]) > 0)
    cat(" column", j, "loads:\n")
    print(signif(N[nz, j], 6))
  }
}
cat("sdr$se_null (outer) columns:", NCOL(sdr$se_null), "\n")
if (NCOL(sdr$se_null)) {
  M <- as.matrix(sdr$se_null); rownames(M) <- on
  for (j in seq_len(ncol(M))) {
    nz <- which(abs(M[, j]) > 0)
    cat(" column", j, ":\n"); print(signif(M[nz, j], 6))
  }
}
# the design: is s(Age)'s fixed column a multiple of Age, or of
# Age and the intercept?
X <- fit$frame$linpreds[[1]]$X
cn <- colnames(X)
cat("mu design columns:", paste(cn, collapse = ", "), "\n")
ia <- match("Age", cn); ifx <- grep("^s[(]Age[)][.]?fx", cn)
if (!is.na(ia) && length(ifx)) {
  fx <- as.numeric(X[, ifx[1]]); ag <- as.numeric(X[, ia])
  r1 <- stats::lm(fx ~ 0 + ag); r2 <- stats::lm(fx ~ ag)
  cat("fx1 ~ 0 + Age: max |resid| =", format(max(abs(resid(r1))), digits = 3),
      "| fx1 ~ 1 + Age: max |resid| =", format(max(abs(resid(r2))), digits = 3),
      "coef", format(coef(r2), digits = 8), "\n")
}
# the warning itself, and the rows' cosines
nd <- fit$data[1:5, ]
w <- character()
invisible(withCallingHandlers(stats::fitted(fit, newdata = nd),
  warning = function(x) { w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning") }))
cat("fitted(newdata) warnings:", length(w), "\n"); for (x in w) cat("  W:", x, "\n")
w <- character()
invisible(withCallingHandlers(stats::fitted(fit, newdata = nd, re_formula = NA),
  warning = function(x) { w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning") }))
cat("fitted(newdata, re_formula = NA) warnings:", length(w), "\n"); for (x in w) cat("  W:", x, "\n")
# finite steps along each outer null column: the objective (inner
# problem re-solved) and the fixed part of mu on the new rows
obj <- fit$obj
p0 <- fit$opt$par
f0 <- obj$fn(p0)
if (NCOL(sdr$se_null)) {
  M <- as.matrix(sdr$se_null)
  u <- fit$par_units %||% rep(1, length(p0))
  for (j in seq_len(ncol(M))) {
    dir <- M[, j] / u
    dir <- dir / sqrt(sum(dir^2))
    for (t in c(1e-3, 1e-2, 1e-1)) {
      f1 <- obj$fn(p0 + t * dir)
      cat(sprintf(" null col %d step %.0e: objective change %.3e\n",
                  j, t, f1 - f0))
    }
    nzb <- which(names(p0) == "beta" & abs(dir) > 0)
    if (length(nzb)) {
      bi <- which(names(p0) == "beta")
      db <- dir[bi]
      Xn <- tryCatch(ns$model_matrix_newdata(fit, nd), error = function(e) NULL)
      cat("  beta loadings of this column:",
          paste0(cn[abs(db) > 0], "=", signif(db[abs(db) > 0], 4),
                 collapse = ", "), "\n")
      # change of X beta on the fit's own first 5 rows along the column
      cat("  d(X beta) on rows 1:5 per unit step:",
          format(as.numeric(as.matrix(X[1:5, , drop = FALSE]) %*% db),
                 digits = 4), "\n")
    }
  }
}
obj$fn(p0)

# The design's exact alias: fx1 = a + b * Age, so d = (Intercept -a,
# Age -b, fx1 1) has X d = 0 on every row, and the likelihood and every
# prediction built from the same basis are constant along it.
if (!is.na(ia) && length(ifx)) {
  cf <- coef(r2)
  bi <- which(names(p0) == "beta")
  d <- numeric(length(p0))
  d[bi[match("(Intercept)", cn)]] <- -cf[1]
  d[bi[ia]] <- -cf[2]
  d[bi[ifx[1]]] <- 1
  cat("max |X d| over rows:", format(max(abs(as.matrix(X) %*% d[bi])),
                                    digits = 3), "\n")
  for (t in c(1e-2, 1e-1, 1)) {
    cat(sprintf(" exact alias step %.0e: objective change %.3e\n", t,
                obj$fn(p0 + t * d) - f0))
  }
  dd <- d / sqrt(sum(d^2))
  M <- as.matrix(sdr$se_null)
  cat("cosine of the exact alias with each lost column:",
      format(abs(colSums(M * dd)) / sqrt(colSums(M^2)), digits = 4), "\n")
  cat("Intercept share of the exact alias (unit norm):",
      format(dd[bi[match("(Intercept)", cn)]], digits = 4), "\n")
  obj$fn(p0)
}
