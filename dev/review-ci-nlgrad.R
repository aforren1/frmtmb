# Reviewer probe for item 2. Two questions the lane's probe does not ask:
#  (a) Do the generators of every family OTHER than the six changed ones
#      give bitwise the same response before and after? (no fits)
#  (b) For the six changed families, how far is each fit from the
#      quantity that raises the ubuntu warning: the maximum absolute
#      gradient (grad_tol 1e-3) and the Newton headroom grad_verdict()
#      computes? A margin of a few is not robust to a platform's rounding.
#
#   Rscript dev/review-ci-nlgrad.R <lib1;lib2;...> <base test> <after test>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[[1L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), "\n")
load_env <- function(path) {
  src <- parse(path)
  keep <- vapply(src, function(e) {
    is.call(e) && identical(e[[1L]], as.name("<-"))
  }, NA)
  env <- new.env(parent = asNamespace("frmtmb"))
  for (e in src[keep]) eval(e, env)
  env
}
eb <- load_env(args[[2L]])
ea <- load_env(args[[3L]])
stopifnot(identical(eb$nl_fams_single, ea$nl_fams_single))
changed <- c("student", "negbinomial", "negbinomial2", "nbinom1",
             "zero_inflated_negbinomial", "hurdle_negbinomial")
same <- vapply(ea$nl_fams_single, function(f) {
  zb <- eb$nl_fam_data(); yb <- eb$nl_fam_response(f, zb)
  rb <- .Random.seed
  za <- ea$nl_fam_data(); ya <- ea$nl_fam_response(f, za)
  identical(yb, ya)
}, NA)
cat("(a) families with bitwise-equal response:", sum(same), "of",
    length(same), "\n")
cat("    differing:", paste(names(same)[!same], collapse = " "), "\n")
cat("    differing but not in the changed six:",
    paste(setdiff(names(same)[!same], changed), collapse = " "), "\n")
cat("    in the changed six but equal:",
    paste(intersect(names(same)[same], changed), collapse = " "), "\n")

fit1 <- function(env, form, f, d) {
  w <- character()
  fit <- withCallingHandlers(
    suppressMessages(frm(form, family = env$nl_fam_of(f), data = d)),
    warning = function(cnd) {
      w <<- c(w, substr(conditionMessage(cnd), 1, 50))
      invokeRestart("muffleWarning")
    })
  v <- frmtmb:::grad_verdict(fit)
  p <- fit$opt$par
  dp <- p[names(p) == "betad"]
  data.frame(ll = as.numeric(logLik(fit)), gmax = v$gmax,
             gpar = v$gmax_par, headroom = v$headroom,
             dpar = paste(sprintf("%.3g", dp), collapse = ","),
             conv = fit$opt$convergence,
             warn = paste(unique(w), collapse = " | "))
}
rows <- list()
for (arm in c("before", "after")) {
  env <- if (arm == "before") eb else ea
  for (f in changed) {
    z <- env$nl_fam_data()
    d <- data.frame(x = z$x, ntr = 10L)
    d$y <- env$nl_fam_response(f, z)
    L <- env$nl_fam_lhs(f)
    fs <- list(
      lin = bf(stats::as.formula(paste(L, "~ x"))),
      int = bf(stats::as.formula(paste(L, "~ 1"))),
      nl = bf(stats::as.formula(paste(L, "~ b0 + b1 * x")), b0 + b1 ~ 1,
              nl = TRUE))
    r <- lapply(fs, function(fm) fit1(env, fm, f, d))
    ratio <- abs(r$nl$ll - r$lin$ll) / (r$lin$ll - r$int$ll)
    for (m in names(r)) {
      rows[[length(rows) + 1L]] <- cbind(arm = arm, family = f, model = m,
                                         r[[m]], gain = r$lin$ll - r$int$ll,
                                         ratio = ratio)
    }
  }
}
out <- do.call(rbind, rows)
op <- options(width = 250)
print(format(out, digits = 4), row.names = FALSE)
cat("\nlargest gmax per arm (grad_tol 1e-3):\n")
print(tapply(out$gmax, out$arm, max))
