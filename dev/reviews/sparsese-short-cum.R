# Short check of lane sparsese, fix round 1, F2: fitted() Est.Error of a
# cumulative REML fit against an analytic delta-method reference written
# here, independent of the package's covariance code.
#   Rscript sparsese-short-cum.R <s30|s4500> <seed> [lib: lane|base]
# s30: dev/sparsese-fix1.R's cumulative design (30 groups, 300 rows,
#      yc ~ x1..x5 + (1 | g)), fitted under frmtmb.joint_root_max 1e6
#      (dense route) and 20 (sparse route) in one process.
# s4500: dev/reviews/sparsese-review-cum4500.R's design (4500 groups of 2
#      rows, 4508 joint parameters), default limit (sparse route).
#
# Reference. Logit cumulative model, eta = X beta + b[g], b ~ N(0, s^2),
# P(Y <= k) = F(c_k - eta), c_1 = tau_1, c_2 = tau_1 + exp(tau_2). Inner
# coordinates (gamma, b), beta = T gamma with MY whitener (QR of the
# scaled design). Outer psi = (tau_1, tau_2, log s). For fixed psi the
# inner mode is found by Newton to machine precision; the REML (Laplace)
# objective is L(psi) = nll(mode) + 0.5 log det H_rr, by block elimination
# (H_bb is diagonal). Vf = inverse of numDeriv's Richardson Hessian of L,
# taken in the coordinates z with psi = psi_hat + U z, U a first-pass
# factor of Vf, so the steps are a hundredth of a standard error.
# J = -W H_r,psi (analytic mixed partials), and per output
#   var = t^2 v'W v + (g_psi - t H_r,psi' W v)' Vf (same),
# with t = dP/deta, v the row's inner load, g_psi = dP/dpsi.
# Where the package lost the group sd (held at its estimate, ~0), the
# reference holds it too: b is fixed at its mode and carries no variance.
a <- commandArgs(trailingOnly = TRUE)
design <- a[1]
seed <- as.integer(a[2])
arm <- if (length(a) >= 3L) a[3] else "lane"
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "design", design, "seed", seed,
    "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
if (design == "s30") {
  set.seed(seed)
  ng <- 30L
  n <- 300L
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(seq_len(ng), length.out = n)),
                  f = factor(rep(c("a", "b"), length.out = n)))
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  d$u <- stats::rnorm(ng, 0, 0.5)[d$g]
  d$e <- stats::rnorm(n, 0, 0.3)
  d$z <- stats::rlogis(n)
  eta <- sin(3 * d$x) + d$u + 0.3 * (d$f == "b")
  d$yc <- factor(cut(0.8 * eta + d$z, c(-Inf, 0, 0.8, Inf),
                     labels = c("lo", "mid", "hi")), ordered = TRUE)
  nd <- d[c(1:6, 151:156), ]
  lims <- c(dense = 1e6, sparse = 20)
} else {
  set.seed(1)
  ng <- 4500L
  n <- 2L * ng
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(seq_len(ng), length.out = n)))
  for (k in 1:5) d[[paste0("x", k)]] <- d$x^k
  eta <- sin(3 * d$x) + stats::rnorm(ng, 0, 0.5)[d$g]
  z <- 0.8 * eta + stats::rlogis(n)
  d$yc <- factor(cut(z, c(-Inf, 0, 0.8, Inf), labels = c("lo", "mid", "hi")),
                 ordered = TRUE)
  nd <- d[1:6, ]
  lims <- c(sparse = NA)
}
fo <- yc ~ x1 + x2 + x3 + x4 + x5 + (1 | g)
on <- function(lim, expr) {
  op <- options(frmtmb.joint_root_max = if (is.na(lim)) NULL else lim)
  on.exit(options(op))
  expr
}
pk <- list()
for (nm in names(lims)) on(lims[[nm]], {
  fit <- q(frm(fo, data = d, family = cumulative(), REML = TRUE))
  sdr <- q(frmtmb:::sdr_of(fit))
  jc <- q(frmtmb:::get_joint_cov(fit))
  cat(sprintf("%s: joint %d  route %s  lost [%s]  rows attr %d\n", nm,
              nrow(sdr$jointPrecision), class(jc$V)[1],
              paste(names(sdr$se_lost), sdr$se_lost, collapse = " "),
              length(attr(sdr$joint_root, "rows"))))
  t0 <- proc.time()[["elapsed"]]
  fn <- q(fitted(fit, newdata = nd))
  fi <- q(fitted(fit))
  cat(sprintf("  fitted() new + in-sample %.2f s\n",
              proc.time()[["elapsed"]] - t0))
  # the drawn set fit_fd_se() perturbs, and fd_joint_root()'s positions
  ds <- frmtmb:::fit_draw_space(fit)
  rn <- jc$names
  pos_d <- unlist(lapply(unique(ds$map$comp), function(cp) which(rn == cp)))
  ov <- frmtmb:::fit_outer_vector(fit, ds$map)
  cat("  drawn (fit_outer_vector):", names(ov), "\n")
  cat("  fd_joint_root pos_d names:", rn[pos_d], " same length:",
      length(pos_d) == length(ov), " root rows:", !is.null(ds$root), "\n")
  b_idx <- seq_along(fit$estimates$b)
  fr <- frmtmb:::fd_joint_root(fit, ds$map, b_idx)
  pk[[nm]] <- list(fit = fit, fn = fn, fi = fi, sdr = sdr, fr = fr,
                   jcV = if (design == "s30") {
                     as.matrix(jc$V[c(pos_d, which(rn == "b")),
                                    c(pos_d, which(rn == "b"))])
                   }, ov = ov)
})
fit <- pk[[1]]$fit
lost <- pk[[1]]$sdr$se_lost
held <- length(lost) > 0L
if (held) stopifnot(identical(names(lost), "theta_1"))

## ---- the reference
X <- function(dd) stats::model.matrix(~ x1 + x2 + x3 + x4 + x5, dd)[, -1]
Xn <- X(d)
p <- ncol(Xn)
sc <- sqrt(colSums(Xn^2) / n)
qx <- qr(Xn / rep(sc, each = n))
stopifnot(qx$rank == p)
Tm <- matrix(0, p, p)
Tm[qx$pivot, ] <- backsolve(qr.R(qx), diag(p)) * sqrt(n) / sc[qx$pivot]
XT <- Xn %*% Tm
gi <- as.integer(d$g)
yi <- as.integer(d$yc)
est <- fit$estimates
psi_hat <- c(est$tau_raw[1], est$tau_raw[2], est$theta[1])
cat("psi_hat (tau_1, tau_2, log s):", format(psi_hat, digits = 10), "\n")
cuts <- function(psi) c(-Inf, psi[1], psi[1] + exp(psi[2]), Inf)
fz <- function(z) ifelse(is.finite(z), stats::dlogis(z), 0)
# stable in the tails: log A from the tail that keeps its digits, and every
# ratio f(.) / A through logs
pieces <- function(eta, psi, y = yi) {
  cc <- cuts(psi)
  u <- cc[y + 1L] - eta
  l <- cc[y] - eta
  lS <- function(z) stats::plogis(z, lower.tail = FALSE, log.p = TRUE)
  lF <- function(z) stats::plogis(z, log.p = TRUE)
  lA <- ifelse(l > 0, lS(l) + log1p(-exp(lS(u) - lS(l))),
               lF(u) + log1p(-exp(lF(l) - lF(u))))
  ru <- ifelse(is.finite(u), exp(stats::dlogis(u, log = TRUE) - lA), 0)
  rl <- ifelse(is.finite(l), exp(stats::dlogis(l, log = TRUE) - lA), 0)
  su <- ifelse(is.finite(u), 1 - 2 * stats::plogis(u), 0)
  sl <- ifelse(is.finite(l), 1 - 2 * stats::plogis(l), 0)
  ge <- ru - rl
  he <- ge^2 - (ru * su - rl * sl)
  du <- ru * su - ge * ru
  dl <- -rl * sl + ge * rl
  # d ge / d c_k, k = 1, 2
  dc <- cbind((y == 1L) * du + (y == 2L) * dl, (y == 2L) * du + (y == 3L) * dl)
  list(nll = -sum(lA), ge = ge, he = he, dc = dc)
}
rs <- function(v) {
  out <- numeric(ng)
  s <- rowsum(v, gi, reorder = TRUE)
  out[as.integer(rownames(s))] <- s
  out
}
# inner mode and Laplace objective at psi; r = (gamma, b)
inner <- function(psi, r0) {
  r <- r0
  s2 <- exp(2 * psi[3])
  fobj <- function(r) {
    b <- if (held) b_fix else r[p + seq_len(ng)]
    pieces(drop(XT %*% r[seq_len(p)]) + b[gi], psi)$nll +
      if (held) 0 else sum(b^2) / (2 * s2)
  }
  for (it in 1:200) {
    gam <- r[seq_len(p)]
    b <- if (held) b_fix else r[p + seq_len(ng)]
    pc <- pieces(drop(XT %*% gam) + b[gi], psi)
    gg <- drop(crossprod(XT, pc$ge))
    Hgg <- crossprod(XT, pc$he * XT)
    if (held) {
      step <- -solve(Hgg, gg)
    } else {
      gb <- rs(pc$ge) + b / s2
      Hgb <- t(rowsum(pc$he * XT, gi, reorder = TRUE))
      hb <- rs(pc$he) + 1 / s2
      S <- Hgg - Hgb %*% (t(Hgb) / hb)
      dg <- -solve(S, gg - Hgb %*% (gb / hb))
      db <- -(gb + drop(crossprod(Hgb, dg))) / hb
      step <- c(dg, db)
    }
    # Newton with backtracking; the objective is convex in r
    f0 <- fobj(r)
    lam <- 1
    repeat {
      rn <- r + lam * step
      f1 <- fobj(rn)
      if ((is.finite(f1) && f1 <= f0 + 1e-12 * abs(f0)) || lam < 1e-10) break
      lam <- lam / 2
    }
    r <- rn
    if (lam == 1 && max(abs(step)) < 1e-13 * (1 + max(abs(r)))) break
  }
  gam <- r[seq_len(p)]
  b <- if (held) b_fix else r[p + seq_len(ng)]
  pc <- pieces(drop(XT %*% gam) + b[gi], psi)
  Hgg <- crossprod(XT, pc$he * XT)
  if (held) {
    return(list(r = r, L = pc$nll + 0.5 * determinant(Hgg)$modulus[[1]],
                pc = pc, S = Hgg, it = it))
  }
  Hgb <- t(rowsum(pc$he * XT, gi, reorder = TRUE))
  hb <- rs(pc$he) + 1 / s2
  S <- Hgg - Hgb %*% (t(Hgb) / hb)
  L <- pc$nll + sum(b^2) / (2 * s2) + ng * psi[3] +
    0.5 * (sum(log(hb)) + determinant(S)$modulus[[1]])
  list(r = r, L = L, pc = pc, S = S, Hgb = Hgb, hb = hb, it = it)
}
gam_pk <- solve(Tm, est$beta)
b_fix <- est$b
r_start <- if (held) gam_pk else c(gam_pk, est$b)
m0 <- inner(psi_hat, r_start)
cat(sprintf("inner Newton iterations %d; mode vs package: beta rel %.2e, b max abs %.2e\n",
            m0$it, max(abs(drop(Tm %*% m0$r[seq_len(p)]) - est$beta)) /
              max(abs(est$beta)),
            if (held) 0 else max(abs(m0$r[p + seq_len(ng)] - est$b))))
fpar <- if (held) 1:2 else 1:3
Lof <- function(ps) {
  full <- psi_hat
  full[fpar] <- ps
  inner(full, m0$r)$L
}
# the coordinates of the second pass: the package's own outer covariance
# where it is usable (it only sets the steps; Vf below does not depend on
# it beyond the differencing error), else a first pass with small steps
Vp <- q(vcov(fit, full = TRUE))
cat("package vcov(full) rows:", rownames(Vp), "\n")
# vcov(full) rows: theta_1, Intercept[1], log(Intercept[2] - Intercept[1])
ord <- if (nrow(Vp) == 3L && rownames(Vp)[1] == "theta_1") c(2L, 3L, 1L)
V0 <- if (!is.null(ord)) Vp[ord, ord][fpar, fpar, drop = FALSE]
U <- tryCatch(t(chol(V0)), error = function(e) NULL)
if (is.null(U)) {
  H1 <- numDeriv::hessian(Lof, psi_hat[fpar],
                          method.args = list(d = 1e-5, eps = 1e-5, r = 4))
  U <- tryCatch(t(chol(solve((H1 + t(H1)) / 2))),
                error = function(e) diag(1 / sqrt(abs(diag(H1))),
                                         length(fpar)))
  cat("U from a first pass\n")
} else cat("U from the package's vcov(full)\n")
Lz <- function(z) Lof(psi_hat[fpar] + drop(U %*% z))
Hz <- numDeriv::hessian(Lz, numeric(length(fpar)),
                        method.args = list(d = 1e-2, eps = 1e-2,
                                           zero.tol = 1, r = 6))
gz <- numDeriv::grad(Lz, numeric(length(fpar)),
                     method.args = list(d = 1e-2, eps = 1e-2,
                                        zero.tol = 1, r = 6))
Vf <- matrix(0, 3, 3)
Vf[fpar, fpar] <- U %*% solve((Hz + t(Hz)) / 2) %*% t(U)
cat("reference: gradient of L at the package's psi, in sd units:",
    format(gz, digits = 3), "\n")
cat("reference Vf (tau_1, tau_2, log s):", format(Vf, digits = 8), "\n")
if (!is.null(V0)) {
  cat(sprintf("reference Vf vs package vcov(full) (scaled max): %.3e\n",
              max(abs(Vf[fpar, fpar] - V0)) /
                max(abs(Vf[fpar, fpar]))))
}

# H_r,psi: d (inner gradient) / d psi at fixed r
pc <- m0$pc
psi <- psi_hat
e2 <- exp(psi[2])
dcp <- cbind(pc$dc[, 1] + pc$dc[, 2], e2 * pc$dc[, 2])
Hrg <- cbind(crossprod(XT, dcp), 0)
if (!held) {
  s2 <- exp(2 * psi[3])
  Hrb <- cbind(apply(dcp, 2, rs), -2 * m0$r[p + seq_len(ng)] / s2)
}
Si <- solve(m0$S)
# per output position: Est.Error of P(Y = k) at rows of design `Xr`
# (natural x1..x5) and group index `gr`
ref_se <- function(Xr, gr, Vf) {
  At <- Xr %*% Tm
  cc <- cuts(psi)
  eta <- drop(At %*% m0$r[seq_len(p)]) +
    if (held) b_fix[gr] else m0$r[p + gr]
  if (held) {
    Xg <- At %*% Si
    qv <- rowSums(At * Xg)
    h <- Xg %*% Hrg
  } else {
    hb <- m0$hb
    Hgb <- m0$Hgb
    Ai <- At - t(Hgb[, gr, drop = FALSE]) / hb[gr]
    Xg <- Ai %*% Si
    cg <- rowSums(t(Hgb[, gr, drop = FALSE]) * Xg)
    qv <- rowSums(At * Xg) + (1 - cg) / hb[gr]
    Pm <- Hgb %*% (Hrb / hb)
    h <- Xg %*% (Hrg - Pm) + Hrb[gr, , drop = FALSE] / hb[gr]
  }
  out <- matrix(NA_real_, nrow(Xr), 3)
  for (k in 1:3) {
    u <- cc[k + 1L] - eta
    l <- cc[k] - eta
    t_ <- -(fz(u) - fz(l))
    dc1 <- (k == 1L) * fz(u) - (k == 2L) * fz(l)
    dc2 <- (k == 2L) * fz(u) - (k == 3L) * fz(l)
    gpsi <- cbind(dc1 + dc2, e2 * dc2, 0)
    w <- gpsi - t_ * h
    out[, k] <- sqrt(t_^2 * qv + rowSums((w %*% Vf) * w))
  }
  out
}
ref_n <- ref_se(X(nd), as.integer(nd$g), Vf)
Vfp <- matrix(0, 3, 3)
if (!is.null(V0)) Vfp[fpar, fpar] <- V0
refp_n <- ref_se(X(nd), as.integer(nd$g), Vfp)
refp_i <- ref_se(Xn, gi, Vfp)
ref_i <- ref_se(Xn, gi, Vf)
rel <- function(x, y) max(abs(x - y) / abs(y))
for (nm in names(pk)) {
  for (w in c("n", "i")) {
    pe <- if (w == "n") pk[[nm]]$fn[, "Est.Error", ] else
      pk[[nm]]$fi[, "Est.Error", ]
    rf <- if (w == "n") ref_n else ref_i
    rp <- if (w == "n") refp_n else refp_i
    cat(sprintf(paste0("%-6s %-9s Est.Error: max |pkg / ref - 1| %.3e;",
                       " with the package's Vf %.3e  (max pkg %.4f,",
                       " max ref %.4f)\n"), nm,
                if (w == "n") "new rows" else "in-sample", rel(pe, rf),
                rel(pe, rp), max(pe), max(rf)))
  }
}
cat("ref new rows Est.Error rows 1-3:\n")
print(round(ref_n[1:3, ], 6))
cat("pkg (", names(pk)[length(pk)], ") new rows rows 1-3:\n")
print(round(pk[[length(pk)]]$fn[1:3, "Est.Error", ], 6))

## ---- F2 decomposition: K K' + Vrest against the reference's natural
## joint covariance over (drawn parameters, b), s30 only
if (design == "s30" && !held) {
  Hgb <- m0$Hgb
  hb <- m0$hb
  H <- rbind(cbind(m0$S + Hgb %*% (t(Hgb) / hb), Hgb),
             cbind(t(Hgb), diag(hb)))
  W <- solve(H)
  Hr <- rbind(Hrg, Hrb)
  J <- -W %*% Hr
  pr <- p + ng
  # natural, in fd_joint_root()'s order: beta, theta (log s), tau_1, tau_2, b
  Mn <- matrix(0, p + 3 + ng, pr + 3)
  Mn[1:p, 1:p] <- Tm
  Mn[p + 1, pr + 3] <- 1
  Mn[p + 2, pr + 1] <- 1
  Mn[p + 3, pr + 2] <- 1
  Mn[p + 3 + seq_len(ng), p + seq_len(ng)] <- diag(ng)
  natref <- function(Vf) {
    Vg <- matrix(0, pr + 3, pr + 3)
    Vg[1:pr, 1:pr] <- W + J %*% Vf %*% t(J)
    Vg[1:pr, pr + 1:3] <- J %*% Vf
    Vg[pr + 1:3, 1:pr] <- t(J %*% Vf)
    Vg[pr + 1:3, pr + 1:3] <- Vf
    Mn %*% Vg %*% t(Mn)
  }
  for (wv in c("own Vf", "package Vf")) {
    Vref <- natref(if (wv == "own Vf") Vf else Vfp)
    sref <- sqrt(diag(Vref))
    sc_err <- function(R) max(abs(R - Vref) / outer(sref, sref))
    cat("-- reference with", wv, "\n")
    for (nm in names(pk)) {
      fr <- pk[[nm]]$fr
      if (is.null(fr)) {
        cat(nm, ": fd_joint_root() NULL\n")
      } else {
        R <- tcrossprod(fr$K) + fr$Vrest
        cat(sprintf(paste0("%-6s fd_joint_root: K %d x %d, Vrest zero off",
                           " the b block: %s; max |KK'+Vrest - Vref| /",
                           " (s_i s_j) %.3e\n"), nm,
                    nrow(fr$K), ncol(fr$K),
                    all(fr$Vrest[1:(p + 3), ] == 0), sc_err(R)))
        if (nm == "sparse") {
          # Vrest = the group effects' covariance given the drawn parameters
          k <- 1:(p + 3)
          Vc <- Vref[-k, -k] - Vref[-k, k] %*% solve(Vref[k, k], Vref[k, -k])
          cat(sprintf(paste0("       Vrest[b, b] vs reference conditional",
                             " cov: %.3e (rel to max)\n"),
                      max(abs(fr$Vrest[-k, -k] - Vc)) / max(abs(Vc))))
        }
      }
      cat(sprintf(paste0("%-6s natural block jc$V[pos, pos] vs reference:",
                         " %.3e (scaled)\n"), nm, sc_err(pk[[nm]]$jcV)))
    }
  }
}
## ---- which outer covariance is right: the reference's Vf against the
## inverse of numDeriv's Richardson Jacobian of the tape's AD gradient (the
## review's method), in the package's order (theta, tau_1, tau_2)
if (!held) {
  obj <- fit$obj
  of <- fit$opt$par
  Hf <- numDeriv::jacobian(function(x) as.numeric(obj$gr(x)), of,
                           method = "Richardson")
  invisible(obj$fn(of))
  Vad <- solve((Hf + t(Hf)) / 2)[c(2, 3, 1), c(2, 3, 1)]
  cat(sprintf(paste0("outer covariance, scaled max |diff|: reference vs",
                     " AD-gradient Richardson %.3e; package vcov(full) vs",
                     " AD-gradient Richardson %.3e\n"),
              max(abs(Vf - Vad) / sqrt(outer(diag(Vad), diag(Vad)))),
              max(abs(Vfp - Vad) / sqrt(outer(diag(Vad), diag(Vad))))))
}
cat(sprintf("peak working set %.0f MB\n",
            ps::ps_memory_info()[["peak_wset"]] / 2^20))
