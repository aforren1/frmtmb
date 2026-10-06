# Item 3: the silent-NaN-SE cases from lane fixes' final review
# (dev/fixes-rev3-flat.R in that worktree), built the same way. Per
# case: the warnings, the SEs, and the outer Hessian's diagnosis:
# optimHess (sdreport's) against the exact AD Hessian obj$he(), rcond,
# the scaled matrix's eigenvalues, and the parameters on its null space.
#   Rscript dev/nanse-item3.R [lib]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
diag_h <- function(f) {
  p <- f$opt$par
  nm <- frmtmb:::outer_par_names(f)
  H1 <- optimHess(p, f$obj$fn, f$obj$gr)
  H2 <- tryCatch(f$obj$he(p), error = function(e) NULL)
  cat("  par:", paste(sprintf("%s=%.4g", nm, p), collapse = " "), "\n")
  cat("  grad:", paste(sprintf("%.2g", f$obj$gr(p)), collapse = " "), "\n")
  cat("  diag optimHess:", paste(signif(diag(H1), 3), collapse = " "), "\n")
  if (!is.null(H2)) {
    cat("  diag he       :", paste(signif(diag(H2), 3), collapse = " "),
        "\n")
  }
  for (H in list(optimHess = H1, he = H2)) {
    if (is.null(H)) next
    if (!all(is.finite(H))) {
      cat("  Hessian not finite in rows:",
          paste(nm[apply(!is.finite(H), 1, any)], collapse = " "), "\n")
      next
    }
    D <- sqrt(abs(diag(H)))
    ok <- D > 0
    S <- (H / outer(D, D))[ok, ok, drop = FALSE]
    e <- eigen(S, TRUE)
    cat("  rcond", signif(rcond(H), 3), "scaled eig",
        paste(signif(e$values, 3), collapse = " "), "\n")
    k <- which.min(abs(e$values))
    v <- e$vectors[, k]
    cat("  min-|eig| dir:", paste(sprintf("%s=%.3f", nm[ok], v)[abs(v) > 0.1],
                                  collapse = " "), "\n")
  }
}
run <- function(lab, expr) {
  cat("\n==", lab, "\n")
  w <- character()
  f <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  se <- withCallingHandlers(sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed)),
                            warning = function(x) {
                              w <<- c(w, paste("[sdr]", conditionMessage(x)))
                              invokeRestart("muffleWarning")
                            })
  cat("  code", f$opt$convergence, "logLik", format(logLik(f), digits = 10),
      "\n  SE:", paste(signif(se, 3), collapse = " "), "\n")
  cat("  warnings:", if (length(w)) paste(substr(w, 1, 150),
                                           collapse = "\n    ") else "none",
      "\n")
  diag_h(f)
  invisible(f)
}
set.seed(41)
x <- runif(300, 0, 1e5)
d1 <- data.frame(x = x, y = 3 * exp(-2e-5 * x) + rnorm(300, 0, 0.05))
f1 <- run("a * exp(b * x), x 0..1e5, b ~ -2e-5",
    frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d1,
        start = list(beta = c(3, -2e-5))))
d1s <- transform(d1, xs = x / 1e5)
run("control: the same with x / 1e5",
    frm(bf(y ~ a * exp(b * xs), a ~ 1, b ~ 1, nl = TRUE), data = d1s,
        start = list(beta = c(3, -2))))
set.seed(44)
d4 <- data.frame(x = runif(200, 0, 3))
d4$y <- 2 * exp(-0.3 * d4$x) + rnorm(200, 0, 0.2)
run("b at its bound lb = 0",
    frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d4,
        start = list(beta = c(2, 0.1)),
        prior = set_prior("normal(0, 5)", nlpar = "b", lb = 0)))
set.seed(45)
d5 <- data.frame(x = rnorm(200))
d5$y <- 1 + 0.5 * d5$x + log(5e-5) + rnorm(200, 0, 0.3)
run("a + b + log(c0), a, b ~ 1 + x, c0 ~ 1 (c0 near 0)",
    frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x, c0 ~ 1, nl = TRUE),
        data = d5, start = list(beta = c(0.5, 0.25, 0.5, 0.25, 5e-5))))
run("same ridge, c0 near 1 (control)",
    frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x, c0 ~ 1, nl = TRUE),
        data = transform(d5, y = y - log(5e-5)),
        start = list(beta = c(0.5, 0.25, 0.5, 0.25, 1))))
d6 <- transform(d5, y = 1 + log(5e-5) * x + rnorm(200, 0, 0.3))
run("identified: a + log(c0) * x, c0 near 0",
    frm(bf(y ~ a + log(c0) * x, a ~ 1, c0 ~ 1, nl = TRUE),
        data = d6, start = list(beta = c(1, 1e-4))))
run("curved ridge: a + log(c0), a ~ 1 + x, c0 near 0",
    frm(bf(y ~ a + log(c0), a ~ 1 + x, c0 ~ 1, nl = TRUE),
        data = d5, start = list(beta = c(1, 0.5, 1e-4))))
