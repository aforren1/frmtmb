# Lane code8 fix round: run a test file with a mutant assigned into the
# installed namespace, to see a new test fail on the defect it pins.
#   Rscript dev/code8-mutant-tests.R <oldprobe|noratio> <test file>
# oldprobe: headroom_checked() as the lane's first build had it (the best
#   of a few straight probes in place of any prediction above grad_tol;
#   review finding 1).
# noratio: se_curvature_real() without its ratio test (review finding 4).
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-code8-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
if (a[1] == "oldprobe") {
  probe <- function(fit, p, g, Hq, M, tol) {
    obj <- fit$obj
    saved <- obj_state_save(obj)
    on.exit({
      obj_state_restore(obj, saved)
      try(obj$gr(p), silent = TRUE)
    }, add = TRUE)
    f <- function(x) {
      tryCatch(as.numeric(se_quiet(obj$fn(x))), error = function(e) NA_real_)
    }
    gq <- drop(crossprod(M, g))
    if (!all(is.finite(gq)) || !all(is.finite(Hq))) return(NA_real_)
    D <- sqrt(abs(diag(Hq)))
    D[!(D > 0)] <- 1
    gs <- gq / D
    if (max(abs(gs)) > sqrt(2 * tol)) return(NA_real_)
    e <- eigen(Hq / outer(D, D), symmetric = TRUE)
    big <- max(abs(e$values))
    cv <- e$values > se_flat_tol * big
    nat <- function(d) drop(M %*% (d / D))
    probes <- list()
    if (any(cv)) {
      V <- e$vectors[, cv, drop = FALSE]
      dn <- -drop(V %*% (drop(crossprod(V, gs)) / e$values[cv]))
      for (t in c(1, 0.5, 0.25, 0.125)) probes[[length(probes) + 1L]] <- t * dn
    }
    dirs <- list(-gs)
    if (any(!cv)) {
      Vo <- e$vectors[, !cv, drop = FALSE]
      dirs[[2L]] <- -drop(Vo %*% crossprod(Vo, gs))
    }
    for (d in dirs) {
      nd <- sqrt(sum(d^2))
      if (!(nd > 0)) next
      for (aa in 10^(-4:1)) probes[[length(probes) + 1L]] <- aa * d / nd
    }
    f0 <- f(p)
    if (!is.finite(f0)) return(NA_real_)
    fx <- vapply(probes, function(d) f(p + nat(d)), 0)
    if (!any(is.finite(fx))) return(NA_real_)
    max(0, f0 - min(fx[is.finite(fx)]))
  }
  chk <- function(fit, p, g, Hq, M, pred, tol) {
    if (isTRUE(pred <= tol)) return(pred)
    r <- tryCatch(probe(fit, p, g, Hq, M, tol), error = function(e) NA_real_)
    if (isTRUE(r <= tol)) r else pred
  }
  environment(probe) <- ns
  environment(chk) <- ns
  assignInNamespace("headroom_checked", chk, "frmtmb")
}
if (a[1] == "noratio") {
  m <- function(fit, p, free, dir_free, lambda) {
    tol <- fit$control$grad_tol %||% 1e-3
    obj <- fit$obj
    saved <- obj_state_save(obj)
    on.exit(obj_state_restore(obj, saved), add = TRUE)
    d <- numeric(length(p))
    d[free] <- dir_free * sqrt(4 * tol / lambda)
    f <- function(x) tryCatch(se_quiet(obj$fn(x)), error = function(e) NA_real_)
    f0 <- f(p)
    full <- f(p + d) + f(p - d) - 2 * f0
    is.finite(full) && full >= 2 * tol
  }
  environment(m) <- ns
  assignInNamespace("se_curvature_real", m, "frmtmb")
}
cat("mutant:", a[1], "\n")
r <- as.data.frame(test_file(a[2], package = "frmtmb",
                             env = test_env("frmtmb"), reporter = "silent"))
print(r[, c("test", "nb", "failed", "error", "skipped")], right = FALSE)
cat("TOTAL failed", sum(r$failed), "error", sum(r$error), "\n")
for (i in which(r$error | r$failed > 0)) {
  res <- attr(r, "tests")
}
x <- test_file(a[2], package = "frmtmb", env = test_env("frmtmb"),
               reporter = "silent")
for (t in x) for (rr in t$results) {
  if (inherits(rr, c("expectation_failure", "expectation_error"))) {
    cat("[", t$test, "]", substr(conditionMessage(rr), 1, 300), "\n")
  }
}
