# Short check of lane sparsese, fix round 1, F1: frm_lp_basis()'s V_root
# recipe and frm_curve() .se on the sparse route against the review's
# independent analytic reference (sparsese-review-indep.R, whose reference
# section is copied below from its line 25 to line 143; one line is changed,
# marked CHANGED: its held-parameter match never matched).
# Note: the copy prints its reference, then the F1 readers below.
#   Rscript sparsese-short-f1int.R <groups> <seed> <group sd> <limit>
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
# CHANGED from the review's copy: se_lost names the parameter "theta_1" and
# opt$par names it "theta", so the copy held nothing; strip the index
keep <- if (length(lost)) !(names(of) %in% sub("_[0-9]+$", "", names(lost))) else
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
cat("package vcov(full):", format(Vfp, digits = 6), " held:", !keep, "
")
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


## ---- F1 readers
rn <- rownames(sdr$jointPrecision)
pk2me <- function(pos) {
  nm <- rn[pos]
  within <- vapply(seq_along(pos), function(i) {
    sum(rn[seq_len(pos[i])] == nm[i])
  }, integer(1))
  ifelse(nm == "beta", within,
         ifelse(nm == "b", nb + within,
                ifelse(nm == "betad", nb + ng + 1L, nb + ng + 2L)))
}
lev <- levels(d$g)
nd <- data.frame(x = seq(1, 2, by = 0.1), g = d$g[1])
for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
nd2 <- nd
nd2$g <- d$g[seq_len(nrow(nd))]
check_lb <- function(lab, lb) {
  cat("--", lab, ": V_root", !is.null(lb$V_root), " V_rest all zero",
      isTRUE(all(lb$V_rest == 0)), " coef_pos", length(lb$coef_pos),
      " extra_var max", max(lb$extra_var), "\n")
  me <- pk2me(lb$coef_pos)
  A <- as.matrix(lb$A)
  B <- matrix(0, nrow(A), p)
  B[, me] <- A
  B <- B %*% M
  rv <- sqrt(ref_var(B))
  rvp <- sqrt(ref_varp(B))
  if (!is.null(lb$V_root)) {
    rec <- sqrt(rowSums((A %*% lb$V_root)^2) +
                  rowSums((A %*% lb$V_rest) * A) + lb$extra_var)
    report(paste(lab, "V_root recipe"), rec, rv, refp = rvp)
    Mc <- M[me, , drop = FALSE]
    for (wv in c("own", "package")) {
      Vn <- Mc %*% (if (wv == "own") Vg else Vgp) %*% t(Mc)
      s <- sqrt(diag(Vn))
      R <- tcrossprod(lb$V_root) + lb$V_rest
      cat(sprintf(paste0("   V_root V_root' + V_rest vs natural reference",
                         " (%s Vf): max |diff| / (s_i s_j) %.3e,",
                         " diag rel %.3e\n"), wv,
                  max(abs(R - Vn) / outer(s, s)),
                  max(abs(diag(R) - s^2) / s^2)))
      ib <- which(rn[lb$coef_pos] == "b")
      if (length(ib)) {
        cat(sprintf(paste0("     b positions: package diag %s, reference",
                           " diag %s; rows' variance min %.3e\n"),
                    paste(format(range(diag(R)[ib]), digits = 4),
                          collapse = ".."),
                    paste(format(range(s[ib]^2), digits = 4),
                          collapse = ".."), min(rec^2)))
      }
    }
  }
  recv <- sqrt(pmax(0, rowSums((A %*% lb$V) * A) + lb$extra_var))
  report(paste(lab, "A V A' (natural V)"), recv, rv, refp = rvp)
  invisible(B)
}
lb0 <- quiet(frm_lp_basis(fit, newdata = nd, re_formula = NA))
B0 <- check_lb("lp_basis noRE", lb0)
lb1 <- quiet(frm_lp_basis(fit, newdata = nd2))
B1 <- check_lb("lp_basis RE", lb1)
if (ng <= 1000L) {
lbi <- quiet(frm_lp_basis(fit))
Bi <- check_lb("lp_basis in-sample", lbi)
rm(Bi)
}
for (w in c("noRE", "RE")) {
  cv <- tryCatch(quiet(if (w == "noRE") {
    frm_curve(fit, newdata = nd, re_formula = NA, simultaneous = FALSE)
  } else frm_curve(fit, newdata = nd2, re_formula = NULL,
                   simultaneous = FALSE)),
  error = function(e) e)
  if (inherits(cv, "error")) {
    cat("frm_curve", w, "ERROR:", substr(conditionMessage(cv), 1, 200), "\n")
  } else {
    B <- if (w == "noRE") B0 else B1
    report(paste("frm_curve .se", w), cv$.se, sqrt(ref_var(B)),
           refp = sqrt(ref_varp(B)))
    cat("   cov_rel_error", format(attr(cv, "check")$cov_rel_error,
                                   digits = 3), "\n")
  }
}
cat(sprintf("peak working set %.0f MB\n",
            ps::ps_memory_info()[["peak_wset"]] / 2^20))
