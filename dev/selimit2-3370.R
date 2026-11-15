# Lane selimit2, backlog 3370: joint_cov_repair() inverts the inner block
# Qrr, which under REML holds beta. Built where it acts: a gaussian REML
# fit on a raw polynomial (beta_whitener() on) whose group sd runs to 0
# (lost, so sdr_rescue() keeps cov_fixed_prop and the repair runs).
#
# The exact covariance of beta (as illse2-p3reml.R defines it) is the
# conditional one at frmtmb's own variance parameters, s2 (X' V^-1 X)^-1,
# plus J V_theta J' over the determined theta. With the sd lost and held,
# and sigma's J zero at a gaussian (beta_hat does not depend on sigma when
# the sd is at 0, and only through the ratio otherwise), the reference is
# the conditional part in 256 bits; J V_theta J' is printed for sigma.
#   Rscript dev/selimit2-3370.R <base|lane> <deg> <seeds>
a <- commandArgs(TRUE)
arm <- a[1]
deg <- as.integer(a[2])
seeds <- eval(parse(text = a[3]))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(Rmpfr)
})
ns <- asNamespace("frmtmb")
cat(sprintf("arm %s | frmtmb %s from %s\n", arm,
            utils::packageVersion("frmtmb"), dirname(find.package("frmtmb"))))
PREC <- 256L
m <- function(x) mpfr(x, PREC)
mM <- function(A) mpfrArray(as.numeric(A), PREC, dim = dim(as.matrix(A)))
num <- function(x) asNumeric(x)
mp_inv <- function(A) {
  if (!inherits(A, "mpfrMatrix")) A <- mM(A)
  p <- nrow(A)
  B <- mM(diag(p))
  for (k in seq_len(p)) {
    piv <- k - 1L + which.max(num(abs(A[k:p, k])))
    if (piv != k) {
      t1 <- A[k, ]; A[k, ] <- A[piv, ]; A[piv, ] <- t1
      t2 <- B[k, ]; B[k, ] <- B[piv, ]; B[piv, ] <- t2
    }
    s <- A[k, k]
    A[k, ] <- A[k, ] / s
    B[k, ] <- B[k, ] / s
    for (i in seq_len(p)[-k]) {
      f <- A[i, k]
      A[i, ] <- A[i, ] - f * A[k, ]
      B[i, ] <- B[i, ] - f * B[k, ]
    }
  }
  B
}
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "s2_err")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
rs <- function(r) {
  if (all(is.na(r))) return("NaN")
  sprintf("%.9f-%.9f (max dev %.1e)", min(r, na.rm = TRUE),
          max(r, na.rm = TRUE), max(abs(r - 1), na.rm = TRUE))
}
xn <- paste0("x", seq_len(deg))
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                              "+ (1 | g)"))
fo0 <- stats::as.formula(paste("~", paste(xn, collapse = " + ")))
for (seed in seeds) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  # no group effect: the REML sd runs to 0
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  X <- stats::model.matrix(fo0, d)
  r <- cap(frm(fo, data = d, REML = TRUE))
  if (inherits(r$v, "s2_err")) {
    cat("seed", seed, "ERROR", r$v, "\n")
    next
  }
  f <- r$v
  sdr <- ns$sdr_of(f)
  lost <- sdr$se_lost
  th <- unname(f$opt$par)
  s2 <- exp(2 * th[1])
  t2 <- exp(2 * th[2])
  # s2 (X' V^-1 X)^-1 with V = I + (t2 / s2) Z Z', in 256 bits
  Xm <- mM(X)
  A <- crossprod(Xm)
  for (ix in split(seq_len(nrow(d)), d$g)) {
    sv <- mpfr2array(colSums(Xm[ix, , drop = FALSE]), dim = c(ncol(X), 1L))
    cg <- m(t2) / (m(s2) + length(ix) * m(t2))
    A <- A - cg * (sv %*% t(sv))
  }
  Vc <- m(s2) * mp_inv(A)
  dg <- sqrt(num(diag(Vc)))
  cat(sprintf("seed %d deg %d | code %d | lost %s | whitened %s | sd %.2e | warn %d\n",
              seed, deg, f$opt$convergence,
              if (length(lost)) paste(names(lost), lost, collapse = ",") else "none",
              !is.null(f$beta_whiten), sqrt(t2), length(r$w)))
  vc <- cap(vcov(f))
  if (!inherits(vc$v, "s2_err")) {
    V <- vc$v
    cat(sprintf("  vcov() SE / exact conditional %s\n",
                rs(sqrt(diag(V))[seq_along(dg)] / dg)))
  }
  nd <- data.frame(x = c(1, 1.25, 1.5, 1.75, 2), g = d$g[1])
  for (k in seq_len(deg)) nd[[xn[k]]] <- nd$x^k
  Am <- mM(stats::model.matrix(fo0, nd))
  vpe <- num(rowSums((Am %*% Vc) * Am))
  fp <- cap(fitted(f, newdata = nd, re_formula = NA))
  if (!inherits(fp$v, "s2_err")) {
    cat(sprintf("  fitted(re_formula = NA) SE / exact conditional %s | warn %d\n",
                rs(fp$v[, "Est.Error"] / sqrt(vpe)), length(fp$w)))
  } else cat("  fitted ERROR", substr(fp$v, 1, 100), "\n")
  # what joint_cov_repair() inverts: Qrr in natural coordinates
  Q <- sdr$jointPrecision
  if (!is.null(Q)) {
    rr <- f$obj$env$random
    ev <- eigen(as.matrix(Q[rr, rr]), symmetric = TRUE, only.values = TRUE)$values
    cat(sprintf("  Qrr: %d x %d, condition number %.2e\n", length(rr),
                length(rr), max(ev) / min(ev)))
  }
}
