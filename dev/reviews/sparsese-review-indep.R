# Review of lane sparsese: the sparse route against a dense joint
# covariance written here, not the package's dense route.
#   Rscript sparsese-review-indep.R <groups> <seed> <group sd> <limit>
# `limit` is options(frmtmb.joint_root_max) for the fit and every reader
# (100 forces the sparse route; "default" leaves 4000).
#
# The reference is analytic. The joint negative log-likelihood of the
# gaussian random-intercept model is
#   ||y - X beta - Z b||^2 / (2 s^2) + ||b||^2 / (2 sd^2) + n log s + q log sd,
# s = exp(betad), sd = exp(theta). In coordinates of MY OWN whitener
# (beta = T gamma, T from a QR of the scaled raw design, not the package's
# beta_whitener()), with C = [X T, Z] and e the residual at the package's
# estimates:
#   H_rr = C'C / s^2 + diag(0, I / sd^2),
#   H_r,betad = 2 C'e / s^2,  H_r,theta = (0, -2 b / sd^2),
#   W = H_rr^-1 (dense chol2inv), J = -W H_rf,
#   Vf = inverse of the REML objective's Hessian over (betad, theta)
#        (numDeriv::jacobian of obj$gr, Richardson); where the fit lost a
#        parameter, the inverse over the others and zero along the lost
#        one (the semantics joint_cov_repair() documents),
#   V_gamma = [W + J Vf J', J Vf; Vf J', Vf].
# A natural load `a` reads (M' a)' V_gamma (M' a), M = T on the beta rows.
# Also an outside check for fitted(re_formula = NA): GLS at the fitted
# variance ratio in the orthogonal basis plus the propagation term.
a <- commandArgs(trailingOnly = TRUE)
ng <- as.integer(a[1])
seed <- as.integer(a[2])
gsd <- as.numeric(a[3])
lim <- a[4]
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "groups", ng, "seed", seed,
    "gsd", gsd, "limit", lim, "\n")
if (lim != "default") options(frmtmb.joint_root_max = as.numeric(lim))
quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (!grepl("a product such as a' V a", m, fixed = TRUE)) {
      cat("   W:", substr(m, 1, 200), "\n")
    }
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
}
set.seed(seed)
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(ng, 0, gsd)[d$g] + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
t0 <- proc.time()[["elapsed"]]
fit <- quiet(frm(fo, data = d, REML = TRUE))
sdr <- quiet(frmtmb:::sdr_of(fit))
cat(sprintf("fit %.2f s\n", proc.time()[["elapsed"]] - t0))
lost <- sdr$se_lost
cat("lost:", paste(names(lost), lost), "\n")
cat("joint_root rows attr:", length(attr(sdr$joint_root, "rows")),
    " route:", class(quiet(frmtmb:::get_joint_cov(fit))$V)[1], "\n")

## ---- the reference
t0 <- proc.time()[["elapsed"]]
X <- function(dd) stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, dd)
Xn <- X(d)
nb <- ncol(Xn)
sc <- sqrt(colSums(Xn^2) / n)
qx <- qr(Xn / rep(sc, each = n))
stopifnot(qx$rank == nb)
Tm <- matrix(0, nb, nb)
Tm[qx$pivot, ] <- backsolve(qr.R(qx), diag(nb)) * sqrt(n) / sc[qx$pivot]
XT <- Xn %*% Tm
of <- fit$opt$par
s2 <- exp(2 * of[["betad"]])
sd2 <- exp(2 * of[["theta"]])
beta_n <- fit$estimates$beta
bhat <- fit$estimates$b
gi <- as.integer(d$g)
e <- d$y - drop(Xn %*% beta_n) - bhat[gi]
p <- nb + ng + 2L
r <- seq_len(nb + ng)
f <- nb + ng + 1:2
jb <- seq_len(nb)
jbb <- nb + seq_len(ng)
Hrr <- matrix(0, nb + ng, nb + ng)
Hrr[jb, jb] <- crossprod(XT) / s2
XtZ <- t(rowsum(XT, gi, reorder = TRUE)) / s2 # nb x ng
Hrr[jb, jbb] <- XtZ
Hrr[jbb, jb] <- t(XtZ)
Hrr[cbind(jbb, jbb)] <- tabulate(gi, ng) / s2 + 1 / sd2
Ce <- c(drop(crossprod(XT, e)), drop(rowsum(e, gi, reorder = TRUE)))
Hrf <- cbind(2 * Ce / s2, c(numeric(nb), -2 * bhat / sd2))
W <- chol2inv(chol(Hrr))
J <- -W %*% Hrf
obj <- fit$obj
Hf <- numDeriv::jacobian(function(x) as.numeric(obj$gr(x)), of,
                         method = "Richardson")
Hf <- (Hf + t(Hf)) / 2
invisible(obj$fn(of))
keep <- if (length(lost)) !(names(of) %in% names(lost)) else
  rep(TRUE, length(of))
Vf <- matrix(0, 2, 2)
Vf[keep, keep] <- solve(Hf[keep, keep, drop = FALSE])
cat("Vf (reference):", format(Vf, digits = 8), "\n")
Vg <- matrix(0, p, p)
Vg[r, r] <- W + J %*% Vf %*% t(J)
Vg[r, f] <- J %*% Vf
Vg[f, r] <- t(Vg[r, f])
Vg[f, f] <- Vf
rm(W)
M <- diag(1, p)
M[jb, jb] <- Tm
cat(sprintf("reference formed %.2f s; inner gradient check max|C'e/s2 - (0,b/sd2)| %.2e\n",
            proc.time()[["elapsed"]] - t0,
            max(abs(Ce / s2 - c(numeric(nb), bhat / sd2)))))
gamma_load <- function(A, pos) {
  B <- matrix(0, nrow(A), p)
  B[, pos] <- A
  B %*% M
}
ref_var <- function(B) rowSums((B %*% Vg) * B)
# the same reference with the package's own Vf (vcov(full = TRUE)), which
# isolates the route from the outer Hessian's own error
Vfp <- unname(quiet(vcov(fit, full = TRUE)))
if (length(lost)) {
  Vfp[!keep, ] <- 0
  Vfp[, !keep] <- 0
}
Vgp <- Vg
JVp <- J %*% Vfp
Vgp[r, r] <- Vg[r, r] - J %*% Vf %*% t(J) + JVp %*% t(J)
Vgp[r, f] <- JVp
Vgp[f, r] <- t(JVp)
Vgp[f, f] <- Vfp
ref_varp <- function(B) rowSums((B %*% Vgp) * B)
rel <- function(x, y) max(abs(x - y) / pmax(abs(x), abs(y)))
relm <- function(x, y) max(abs(x - y)) / max(abs(y))
report <- function(nm, pkg, ref, how = rel, refp = NULL) {
  cat(sprintf("%-26s n=%5d rel %.3e%s\n", nm, length(ref), how(pkg, ref),
              if (is.null(refp)) "" else
                sprintf("  (package Vf: %.3e)", how(pkg, refp))))
}

## ---- readers
lev <- levels(d$g)
nd <- data.frame(x = seq(1, 2, by = 0.1), g = d$g[1])
for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
nd2 <- nd
nd2$g <- d$g[seq_len(nrow(nd))]
A0 <- X(nd)
t0 <- proc.time()[["elapsed"]]
se <- quiet(fitted(fit, newdata = nd, re_formula = NA))
t1 <- proc.time()[["elapsed"]] - t0
B0 <- gamma_load(A0, jb)
v0 <- ref_var(B0)
report("fitted noRE", se[, "Est.Error"], sqrt(v0), refp = sqrt(ref_varp(B0)))
cat(sprintf("   (fitted noRE %.2f s; mean check %.1e)\n", t1,
            rel(se[, "Estimate"], drop(A0 %*% beta_n))))
# outside check of the conditional part: GLS, orthogonal basis
tt <- sd2 / s2
pb <- stats::poly(d$x, 7)
Xo <- cbind(1, pb)
Xt <- Xo
for (ix in split(seq_len(n), d$g)) {
  m <- length(ix)
  aa <- (1 - 1 / sqrt(1 + m * tt)) / m
  Xt[ix, ] <- Xo[ix, , drop = FALSE] -
    aa * rep(colSums(Xo[ix, , drop = FALSE]), each = m)
}
Ri <- backsolve(qr.R(qr(Xt)), diag(ncol(Xt)))
Ao <- cbind(1, stats::predict(pb, nd$x))
cond <- s2 * rowSums((Ao %*% Ri)^2)
BJ <- B0[, r] %*% J
prop_ref <- rowSums((BJ %*% Vf) * BJ)
report("fitted noRE vs GLS+prop", se[, "Est.Error"], sqrt(cond + prop_ref))
report("reference vs GLS+prop", sqrt(v0), sqrt(cond + prop_ref))
cat("   propagation share of variance:",
    format(range(prop_ref / v0), digits = 3), "\n")
A1 <- cbind(X(nd2), 1)
B1 <- t(vapply(seq_len(nrow(nd2)), function(i) {
  drop(gamma_load(A1[i, , drop = FALSE],
                  c(jb, jbb[match(as.character(nd2$g[i]), lev)])))
}, numeric(p)))
se1 <- quiet(fitted(fit, newdata = nd2))
report("fitted RE", se1[, "Est.Error"], sqrt(ref_var(B1)),
       refp = sqrt(ref_varp(B1)))
cat(sprintf("   (mean check %.1e)\n", rel(se1[, "Estimate"],
    drop(X(nd2) %*% beta_n) + bhat[match(as.character(nd2$g), lev)])))
ln <- quiet(frm_linpred(fit, newdata = nd2, se.fit = TRUE))
report("frm_linpred se.fit", ln$se.fit, sqrt(ref_var(B1)),
       refp = sqrt(ref_varp(B1)))
t0 <- proc.time()[["elapsed"]]
fi <- quiet(fitted(fit))
t1 <- proc.time()[["elapsed"]] - t0
Bi <- matrix(0, n, p)
Bi[, jb] <- XT
Bi[cbind(seq_len(n), jbb[gi])] <- 1
report("fitted in-sample", fi[, "Est.Error"], sqrt(ref_var(Bi)),
       refp = sqrt(ref_varp(Bi)))
cat(sprintf("   (in-sample %.2f s; mean check %.1e)\n", t1,
            rel(fi[, "Estimate"], drop(Xn %*% beta_n) + bhat[gi])))
rm(Bi)
Vn <- Tm %*% Vg[jb, jb] %*% t(Tm)
Vnp <- Tm %*% Vgp[jb, jb] %*% t(Tm)
V <- quiet(vcov(fit))
report("vcov (max entry)", unname(V), Vn, relm, refp = Vnp)
report("vcov diag", sqrt(diag(V)), sqrt(diag(Vn)), refp = sqrt(diag(Vnp)))
Vfull <- quiet(vcov(fit, full = TRUE))
cat("vcov(full) dimnames:", rownames(Vfull), "\n")
report("vcov(full)", unname(Vfull), Vf, relm)
Cm <- rbind(c(0, 1, 1, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0, 0, 1),
            c(1, 0, 0, 0, 0, 0, 0, 0))
h <- quiet(hypothesis(fit, c("x1 + x2 = 0", "x7 = 0", "Intercept = 0")))
report("hypothesis", h$hypothesis$Est.Error,
       sqrt(ref_var(gamma_load(Cm, jb))),
       refp = sqrt(ref_varp(gamma_load(Cm, jb))))
dcr <- sdr$diag.cov.random
pr <- names(sdr$par.random)
if (length(dcr) == length(pr)) {
  report("diag.cov.random beta", dcr[pr == "beta"], diag(Vn),
         refp = diag(Vnp))
  report("diag.cov.random b", dcr[pr == "b"], diag(Vg)[jbb],
         refp = diag(Vgp)[jbb])
}
re <- tryCatch(quiet(ranef(fit)), error = function(e) e)
if (!inherits(re, "error")) {
  re2 <- tryCatch(quiet(ranef(fit, summary = TRUE)), error = function(e) e)
  rg <- if (inherits(re2, "error")) re$g else re2$g
  cat("   ranef class/dims:", class(rg), dim(rg), "attrs:",
      names(attributes(rg)), "\n")
  cv <- attr(rg, "condVar") %||% attr(rg, "postVar")
  ese <- if (!is.null(cv)) sqrt(as.numeric(cv)) else
    if (length(dim(rg)) == 3L && dim(rg)[2] >= 2L) rg[, 2L, 1L] else NULL
  if (!is.null(ese)) {
    report("ranef Est.Error", unname(ese), sqrt(diag(Vg)[jbb]),
           refp = sqrt(diag(Vgp)[jbb]))
  }
} else cat("ranef error:", conditionMessage(re), "\n")
cat(sprintf("peak working set %.0f MB\n",
            ps::ps_memory_info()[["peak_wset"]] / 2^20))
