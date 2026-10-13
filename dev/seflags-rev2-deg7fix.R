# Reviewer of lane seflags, claim 4: can frmtmb's own objective give the
# degree-7 SEs? The finite-difference outer Hessian cannot resolve its
# smallest directions (noise 1e-12 against eigenvalues 1e-11 to 1e-14 of
# the largest), but the OBJECTIVE resolves them: re-measure the Hessian
# restricted to the m smallest unit-diagonal eigenvectors by symmetric
# second differences of obj$fn at steps sized to move it by about 1
# (the step test's own kind of evaluation), splice that block into the
# eigen-decomposition, invert, and compare with lmer()'s SEs.
#   Rscript dev/seflags-rev-deg7fix.R <lane|base> [seed]
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev2-common.R")
a <- commandArgs(trailingOnly = TRUE)
sd0 <- if (length(a) > 1) as.integer(a[2]) else 7L
deg <- 7L
set.seed(sd0)
d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
  stats::rnorm(200, 0, 0.2)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
rhs <- paste0("x", 1:deg, collapse = " + ")
fo <- stats::as.formula(paste("y ~", rhs, "+ (1 | g)"))
cn <- c("(Intercept)", paste0("x", 1:deg))
m2 <- suppressWarnings(lme4::lmer(fo, data = d, REML = FALSE,
  control = lme4::lmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun = 1e5))))
ref <- sqrt(diag(as.matrix(vcov(m2))))
f0 <- suppressWarnings(frm(fo, data = d))
on <- ns$outer_par_names(f0)
p <- f0$opt$par
q <- p
q[match(cn, on)] <- lme4::fixef(m2)
q[on == "theta_1"] <- log(attr(lme4::VarCorr(m2)$g, "stddev")[[1]])
q[grepl("^sigma", on)] <- log(stats::sigma(m2))
fn <- function(x) f0$obj$fn(x)
for (pt in list(list("frmtmb's optimum", p), list("lmer's optimum", q))) {
  x0 <- pt[[2]]
  h <- ns$fd_hessian(x0, f0$obj$gr)
  H <- h$H
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), symmetric = TRUE)
  np <- length(x0)
  c2 <- function(dir) {
    # scale so the quadratic model moves fn by about 1, then measure
    s <- 1
    f00 <- fn(x0)
    for (it in 1:60) {
      v <- fn(x0 + s * dir) + fn(x0 - s * dir) - 2 * f00
      if (is.finite(v) && v > 0.5) break
      s <- s * 4
    }
    v / s^2
  }
  for (m in 2:3) {
    idx <- (np - m + 1):np
    V <- e$vectors[, idx, drop = FALSE]
    # directions in natural units
    U <- V / D
    B <- matrix(0, m, m)
    for (i in seq_len(m)) B[i, i] <- c2(U[, i])
    for (i in seq_len(m)) for (j in seq_len(m)) if (i < j) {
      cij <- c2(U[, i] + U[, j])
      B[i, j] <- B[j, i] <- (cij - B[i, i] - B[j, j]) / 2
    }
    S <- e$vectors[, -idx, drop = FALSE] %*%
      (e$values[-idx] * t(e$vectors[, -idx, drop = FALSE])) + V %*% B %*% t(V)
    Vc <- solve(S) / outer(D, D)
    se <- sqrt(diag(Vc))[match(cn, on)]
    cat(sprintf("%s, %d smallest re-measured: ev %s -> %s | SE/lmer %s\n",
                pt[[1]], m, paste(signif(e$values[idx], 3), collapse = " "),
                paste(signif(eigen(B, only.values = TRUE)$values, 3),
                      collapse = " "),
                paste(format(se / ref, digits = 4), collapse = " ")))
  }
}
fn(p)
