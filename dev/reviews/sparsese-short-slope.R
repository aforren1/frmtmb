# Short check of lane sparsese, fix round 1, F1 on a random-slope design:
# frm_lp_basis()'s V_root recipe, V_root V_root' + V_rest, and frm_curve()
# .se against an analytic reference written here.
#   Rscript sparsese-short-slope.R <seed>
# Design: sparsese-review-designs.R "slope" (30 groups, 300 rows,
# y ~ x1..x5 + (1 + x1 | g), gaussian, REML), fitted under
# frmtmb.joint_root_max 1e6 (dense route) and 20 (sparse route).
#
# Reference. Joint nll in (gamma, b) for fixed psi = (log s, log sd_0,
# log sd_1, z), rho = z / sqrt(1 + z^2) (the package's coordinates, so its
# vcov(full) can stand in for Vf in a second column; the reference does not
# depend on the choice):
#   ||y - XT gamma - Z b||^2 / (2 s^2) + n log s
#     + sum_j b_j' Sigma^-1 b_j / 2 + (ng / 2) log |Sigma|,
# beta = T gamma with MY whitener (QR of the scaled design). H_rr is
# constant in r, the mode is one solve, the REML objective is exact:
# L(psi) = nll(mode) + 0.5 log det H_rr. Vf = inverse of numDeriv's
# Richardson Hessian of L in coordinates scaled by a standard error;
# J = -W H_r,psi with analytic d Sigma^-1 / d psi;
# V_gamma = [W + J Vf J', J Vf; Vf J', Vf]; natural = M gamma.
a <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(a[1])
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "seed", seed, "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(seed)
ng <- 30L
n <- 300L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)),
                h = factor(sample(rep(seq_len(12), length.out = n))),
                f = factor(rep(c("a", "b"), each = 1, length.out = n)))
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
ug <- stats::rnorm(ng, 0, 0.5)
sg <- stats::rnorm(ng, 0, 0.4)
uh <- stats::rnorm(12, 0, 0.3)
eta <- sin(3 * d$x) + ug[d$g] + 0.3 * (d$f == "b")
d$y <- eta + sg[d$g] * d$x + stats::rnorm(n, 0, 0.3)
fo <- y ~ x1 + x2 + x3 + x4 + x5 + (1 + x1 | g)
lev <- levels(d$g)
nd <- data.frame(x = seq(1, 2, by = 0.1), g = d$g[1])
for (k in 1:5) nd[[paste0("x", k)]] <- nd$x^k
nd2 <- nd
nd2$g <- d$g[seq_len(nrow(nd))]
pk <- list()
for (lim in c(dense = 1e6, sparse = 20)) {
  options(frmtmb.joint_root_max = lim)
  fit <- q(frm(fo, data = d, REML = TRUE))
  sdr <- q(frmtmb:::sdr_of(fit))
  jc <- q(frmtmb:::get_joint_cov(fit))
  nm <- if (lim > 1e5) "dense" else "sparse"
  cat(sprintf("%s: joint %d  whitened %s  route %s  lost [%s]\n", nm,
              nrow(sdr$jointPrecision), !is.null(fit$beta_whiten),
              class(jc$V)[1],
              paste(names(sdr$se_lost), collapse = " ")))
  out <- list(fit = fit, rn = rownames(sdr$jointPrecision),
              lb0 = q(frm_lp_basis(fit, newdata = nd, re_formula = NA)),
              lb1 = q(frm_lp_basis(fit, newdata = nd2)),
              lbi = q(frm_lp_basis(fit)),
              lp1 = q(frm_linpred(fit, newdata = nd2, se.fit = TRUE))$se.fit,
              Vp = q(vcov(fit, full = TRUE)))
  for (w in c("cv0", "cv1")) {
    out[[w]] <- tryCatch(q(frm_curve(fit, newdata = if (w == "cv0") nd else
      nd2, re_formula = if (w == "cv0") NA else NULL,
      simultaneous = FALSE)), error = function(e) e)
  }
  pk[[nm]] <- out
}
options(frmtmb.joint_root_max = NULL)
fit <- pk$dense$fit
cat("estimates identical across routes:",
    identical(fit$estimates, pk$sparse$fit$estimates), "\n")

## ---- the reference
X <- function(dd) stats::model.matrix(~ x1 + x2 + x3 + x4 + x5, dd)
Xn <- X(d)
nb <- ncol(Xn)
sc <- sqrt(colSums(Xn^2) / n)
qx <- qr(Xn / rep(sc, each = n))
stopifnot(qx$rank == nb)
Tm <- matrix(0, nb, nb)
Tm[qx$pivot, ] <- backsolve(qr.R(qx), diag(nb)) * sqrt(n) / sc[qx$pivot]
XT <- Xn %*% Tm
gi <- as.integer(d$g)
nq <- 2L * ng
Zfor <- function(dd) {
  Z <- matrix(0, nrow(dd), nq)
  j <- as.integer(dd$g)
  Z[cbind(seq_len(nrow(dd)), 2L * j - 1L)] <- 1
  Z[cbind(seq_len(nrow(dd)), 2L * j)] <- dd$x1
  Z
}
C <- cbind(XT, Zfor(d))
CtC <- crossprod(C)
Cty <- drop(crossprod(C, d$y))
Sig <- function(ps) {
  s0 <- exp(ps[2]); s1 <- exp(ps[3]); rho <- ps[4] / sqrt(1 + ps[4]^2)
  matrix(c(s0^2, rho * s0 * s1, rho * s0 * s1, s1^2), 2, 2)
}
dSig <- function(ps) {
  s0 <- exp(ps[2]); s1 <- exp(ps[3]); rho <- ps[4] / sqrt(1 + ps[4]^2)
  dr <- (1 + ps[4]^2)^(-1.5)
  list(matrix(c(2 * s0^2, rho * s0 * s1, rho * s0 * s1, 0), 2, 2),
       matrix(c(0, rho * s0 * s1, rho * s0 * s1, 2 * s1^2), 2, 2),
       matrix(c(0, dr * s0 * s1, dr * s0 * s1, 0), 2, 2))
}
pr <- nb + nq
Hof <- function(ps) {
  H <- CtC / exp(2 * ps[1])
  Si <- solve(Sig(ps))
  for (j in seq_len(ng)) {
    k <- nb + 2L * j - 1:0
    H[k, k] <- H[k, k] + Si
  }
  H
}
Lof <- function(ps) {
  H <- Hof(ps)
  R <- chol(H)
  r <- backsolve(R, forwardsolve(t(R), Cty / exp(2 * ps[1])))
  e <- d$y - drop(C %*% r)
  bm <- matrix(r[nb + seq_len(nq)], 2)
  Si <- solve(Sig(ps))
  sum(e^2) / (2 * exp(2 * ps[1])) + n * ps[1] +
    sum(bm * (Si %*% bm)) / 2 + ng / 2 * determinant(Sig(ps))$modulus[[1]] +
    sum(log(diag(R)))
}
psi <- unname(fit$opt$par)
cat("psi_hat (log s, log sd0, log sd1, z):", format(psi, digits = 10), "\n")
H <- Hof(psi)
W <- chol2inv(chol(H))
r <- drop(W %*% Cty) / exp(2 * psi[1])
cat(sprintf("mode vs package: beta rel %.2e\n",
            max(abs(drop(Tm %*% r[seq_len(nb)]) - fit$estimates$beta)) /
              max(abs(fit$estimates$beta))))
# the package's b order, matched to mine by value
bpk <- fit$estimates$b
bme <- r[nb + seq_len(nq)]
bmap <- vapply(bpk, function(v) which.min(abs(bme - v)), integer(1))
cat(sprintf("b matched: one-to-one %s, max |diff| %.2e\n",
            !anyDuplicated(bmap), max(abs(bme[bmap] - bpk))))
cat("vcov(full) rows:", rownames(pk$dense$Vp), "; opt par:", names(fit$opt$par), "
")
Vp <- unname(pk$dense$Vp)
# a repaired fit propagates cov_fixed_prop, not what vcov(full) prints
Vprop <- q(frmtmb:::sdr_of(fit))$cov_fixed_prop
if (!is.null(Vprop)) {
  cat("vcov(full) diag:", format(diag(Vp), digits = 4), "; cov_fixed_prop diag:",
      format(diag(Vprop), digits = 4), "
")
  Vp <- unname(as.matrix(Vprop))
}
# a lost parameter is held at its estimate (joint_cov_repair()'s
# semantics): no variance, and the Hessian over the others
lostn <- names(q(frmtmb:::sdr_of(fit))$se_lost)
hold <- 1L + as.integer(sub("^theta_", "", lostn[grepl("^theta_", lostn)]))
fp <- setdiff(1:4, hold)
cat("held positions of psi:", hold, "\n")
Vp[hold, ] <- 0
Vp[, hold] <- 0
U <- t(chol(Vp[fp, fp]))
Lz <- function(z) {
  ps <- psi
  ps[fp] <- ps[fp] + drop(U %*% z)
  Lof(ps)
}
Hz <- numDeriv::hessian(Lz, numeric(length(fp)),
                        method.args = list(d = 1e-2, eps = 1e-2,
                                           zero.tol = 1, r = 6))
gz <- numDeriv::grad(Lz, numeric(length(fp)),
                     method.args = list(d = 1e-2, eps = 1e-2, zero.tol = 1,
                                        r = 6))
cat("gradient of L at the package's psi, in sd units:", format(gz, digits = 3),
    "\n")
Vf <- matrix(0, 4, 4)
Vf[fp, fp] <- U %*% solve((Hz + t(Hz)) / 2) %*% t(U)
cat(sprintf("Vf: reference vs package vcov(full), scaled max %.3e\n",
            max(abs(Vf - Vp)[fp, fp] / sqrt(outer(diag(Vf), diag(Vf)))[fp, fp])))
e <- d$y - drop(C %*% r)
Hr <- matrix(0, pr, 4)
Hr[, 1] <- 2 * drop(crossprod(C, e)) / exp(2 * psi[1])
Si <- solve(Sig(psi))
dS <- dSig(psi)
bm <- matrix(r[nb + seq_len(nq)], 2)
for (k in 1:3) {
  dSi <- -Si %*% dS[[k]] %*% Si
  Hr[nb + seq_len(nq), k + 1L] <- as.vector(dSi %*% bm)
}
J <- -W %*% Hr
Vgam <- function(Vf) {
  V <- matrix(0, pr + 4, pr + 4)
  V[1:pr, 1:pr] <- W + J %*% Vf %*% t(J)
  V[1:pr, pr + 1:4] <- J %*% Vf
  V[pr + 1:4, 1:pr] <- t(J %*% Vf)
  V[pr + 1:4, pr + 1:4] <- Vf
  V
}
Vg <- Vgam(Vf)
Vgp <- Vgam(Vp)
M <- diag(1, pr + 4)
M[1:nb, 1:nb] <- Tm
# package joint position -> reference position
pk2me <- function(pos, rn) {
  nm <- rn[pos]
  within <- vapply(seq_along(pos), function(i) {
    sum(rn[seq_len(pos[i])] == nm[i])
  }, integer(1))
  out <- integer(length(pos))
  out[nm == "beta"] <- within[nm == "beta"]
  out[nm == "b"] <- nb + bmap[within[nm == "b"]]
  out[nm == "betad"] <- pr + 1L
  out[nm == "theta"] <- pr + 1L + within[nm == "theta"]
  out
}
rel <- function(x, y) max(abs(x - y) / abs(y))
for (nm in names(pk)) {
  o <- pk[[nm]]
  cat("==", nm, "\n")
  for (w in c("lb0", "lb1", "lbi")) {
    lb <- o[[w]]
    me <- pk2me(lb$coef_pos, o$rn)
    A <- as.matrix(lb$A)
    B <- matrix(0, nrow(A), pr + 4)
    B[, me] <- A
    B <- B %*% M
    rv <- sqrt(rowSums((B %*% Vg) * B))
    rvp <- sqrt(rowSums((B %*% Vgp) * B))
    cat(sprintf("  %s: V_root %s, V_rest zero %s, %d rows, %d positions\n",
                w, !is.null(lb$V_root), isTRUE(all(lb$V_rest == 0)),
                nrow(A), length(me)))
    if (!is.null(lb$V_root)) {
      rec <- sqrt(rowSums((A %*% lb$V_root)^2) +
                    rowSums((A %*% lb$V_rest) * A) + lb$extra_var)
      cat(sprintf(paste0("    V_root recipe vs reference: %.3e;",
                         " with the package's Vf %.3e\n"),
                  rel(rec, rv), rel(rec, rvp)))
      Mc <- M[me, , drop = FALSE]
      R <- tcrossprod(lb$V_root) + lb$V_rest
      for (wv in c("own", "package")) {
        Vn <- Mc %*% (if (wv == "own") Vg else Vgp) %*% t(Mc)
        s <- sqrt(diag(Vn))
        cat(sprintf(paste0("    V_root V_root' + V_rest vs natural",
                           " reference (%s Vf): %.3e (scaled max)\n"),
                    wv, max(abs(R - Vn) / outer(s, s))))
      }
    }
    recv <- sqrt(pmax(0, rowSums((A %*% lb$V) * A) + lb$extra_var))
    cat(sprintf("    A V A' (natural V) vs reference (package Vf): %.3e\n",
                rel(recv, rvp)))
    if (w == "lb1") {
      cat(sprintf("    frm_linpred se.fit (RE) vs reference (package Vf): %.3e\n",
                  rel(o$lp1, rvp)))
    }
    for (cw in c("cv0", "cv1")) {
      if ((cw == "cv0") != (w == "lb0") || w == "lbi") next
      cv <- o[[cw]]
      if (inherits(cv, "error")) {
        cat("    frm_curve ERROR:", substr(conditionMessage(cv), 1, 200), "\n")
      } else {
        cat(sprintf(paste0("    frm_curve .se vs reference: %.3e; with the",
                           " package's Vf %.3e; cov_rel_error %.2e\n"),
                    rel(cv$.se, rv), rel(cv$.se, rvp),
                    attr(cv, "check")$cov_rel_error))
      }
    }
  }
}
cat(sprintf("peak working set %.0f MB\n",
            ps::ps_memory_info()[["peak_wset"]] / 2^20))
