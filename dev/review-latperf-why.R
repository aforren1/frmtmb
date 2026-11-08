# Review of lane latperf: why a fit leaves the sparse route, and what the
# repair branch changes (the norepair mutant survives every test file).
# Data: dev/review-latperf-algebra.R's mk(), seed 20261106.
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat(sprintf("lib: %s frmtmb %s\n", dirname(find.package("frmtmb")),
            format(packageVersion("frmtmb"))))
ns <- asNamespace("frmtmb")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-latperf/dev/review-latperf-algebra.R")
eval(parse(text = src[grep("^mk <- function", src):
                        (grep("^designs <- list", src) - 1L)]))
d <- mk()
q <- function(e) suppressWarnings(suppressMessages(e))

step_through <- function(fit) {
  Q <- frmtmb:::joint_precision(fit)
  sdr <- frmtmb:::sdr_of(fit)
  cat("  cov_fixed_root", !is.null(sdr$cov_fixed_root), "joint_root",
      !is.null(sdr$joint_root), "prop", !is.null(sdr$cov_fixed_prop), "\n")
  p <- nrow(Q)
  r <- fit$obj$env$random
  f <- setdiff(seq_len(p), r)
  Q <- methods::as(methods::as(Q, "CsparseMatrix"), "generalMatrix")
  Qrr <- Matrix::forceSymmetric(Q[r, r, drop = FALSE])
  L <- tryCatch(Matrix::Cholesky(Qrr, perm = TRUE, super = TRUE, LDL = FALSE),
                error = function(e) conditionMessage(e))
  if (is.character(L)) {
    cat("  Cholesky(Q_rr) fails:", L, "\n")
    ev <- eigen(as.matrix(Qrr), only.values = TRUE)$values
    cat(sprintf("  Q_rr eigen min %.3g max %.3g\n", min(ev), max(ev)))
    return(invisible())
  }
  J <- -as.matrix(Matrix::solve(L, Q[r, f, drop = FALSE]))
  cat("  J finite", all(is.finite(J)), "\n")
  S <- as.matrix(Q[f, f]) + as.matrix(Matrix::crossprod(Q[r, f], J))
  ok <- tryCatch({solve(S); "ok"}, error = function(e) conditionMessage(e))
  cat("  base solve(S):", ok, " rcond", format(rcond(S), digits = 3), "\n")
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("  S eigen min %.3g max %.3g\n", min(ev), max(ev)))
  dn <- tryCatch({Matrix::solve(Q); "ok"}, error = function(e) conditionMessage(e))
  cat("  Matrix::solve(Q) (the dense route):", substr(dn, 1, 100), "\n")
}

cat("\n== unstr, the algebra script's design\n")
fit <- q(frm(cnt0 ~ x + (1 | id) + unstr(t, g), family = poisson(), data = d))
cat("  code", fit$opt$convergence, "\n")
step_through(fit)
v <- q(fitted(fit))
cat("  fitted Est.Error NA:", sum(is.na(v[, "Est.Error"])), "of", nrow(v), "\n")

cat("\n== what the repair branch changes: boundary and ar2 designs\n")
for (spec in list(
  list("boundary", cnt0 ~ x + (1 | id) + ar(t, g, cov = TRUE)),
  list("ar2", cnt0 ~ x + (1 | id) + ar(t, g, p = 2)))) {
  fit <- q(frm(spec[[2]], family = poisson(), data = d))
  sdr <- frmtmb:::sdr_of(fit)
  a <- q(fitted(fit))[, "Est.Error"]
  old <- get("jcs_joint_cov", envir = ns)
  txt <- deparse(old)
  txt <- sub("repair <- !is.null(Vp) &&", "repair <- FALSE && !is.null(Vp) &&",
             txt, fixed = TRUE)
  mut <- eval(parse(text = txt))
  environment(mut) <- ns
  unlockBinding("jcs_joint_cov", ns)
  assign("jcs_joint_cov", mut, envir = ns)
  f2 <- fit
  f2$cache <- list2env(mget(setdiff(ls(fit$cache), "Vjoint"),
                            envir = fit$cache))
  b <- q(fitted(f2))[, "Est.Error"]
  sp <- frmtmb:::jc_is_sparse(frmtmb:::get_joint_cov(f2))
  assign("jcs_joint_cov", old, envir = ns)
  ok <- is.finite(a) & is.finite(b)
  cat(sprintf(paste("  %s: lost %s; without the repair: sparse %s, rows",
                    "finite %d vs %d, max |rel diff| %.3g\n"),
              spec[[1]], paste(names(sdr$se_lost), collapse = ","), sp,
              sum(is.finite(b)), sum(is.finite(a)),
              if (any(ok)) max(abs(b[ok] / a[ok] - 1)) else NA))
}
