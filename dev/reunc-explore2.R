# Lane wt-reunc: the joint precision's row layout under each estimation
# mode, and the sparse Cholesky solve systems a conditional draw needs.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
G <- 8; m <- 5
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[d$g] + rnorm(G * m)
for (mode in c("ML", "REML", "profile")) {
  fit <- switch(mode,
    ML = frm(bf(y ~ x + (1 | g)), data = d),
    REML = frm(bf(y ~ x + (1 | g)), data = d, REML = TRUE),
    profile = frm(bf(y ~ x + (1 | g)), data = d,
                  control = frmtmb_control(profile = TRUE)))
  Q <- frmtmb:::autoscale_sdreport(fit, jp = TRUE)$jointPrecision
  cat(mode, "Q rows:", paste(names(table(rownames(Q))), table(rownames(Q))),
      " class", class(Q), "\n")
  ds <- frmtmb:::fit_draw_space(fit)
  cat("  draw space:", ds$map$comp, "\n")
  cat("  opt par names:", names(fit$obj$par), "\n")
  cat("  random:", fit$obj$env$random, "\n")
}
Q <- as(Matrix::forceSymmetric(Q), "CsparseMatrix")
ch <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)
z <- rnorm(nrow(Q))
x1 <- Matrix::solve(ch, Matrix::solve(ch, z, system = "Lt"), system = "Pt")
print(class(x1))
# check: covariance of P' L^-T z is Q^-1
Z <- matrix(rnorm(nrow(Q) * 20000), nrow(Q))
X <- as.matrix(Matrix::solve(ch, Matrix::solve(ch, Z, system = "Lt"),
                             system = "Pt"))
Vh <- tcrossprod(X) / ncol(X)
V <- as.matrix(Matrix::solve(Q))
cat("max rel diff diag:", max(abs(diag(Vh) / diag(V) - 1)), "\n")
