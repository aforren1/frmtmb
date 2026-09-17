## Reviewer, claim 2: VarCorr(fit) Est.Error against an independent delta
## method. Independent parts: the Jacobian is numDeriv's Richardson
## extrapolation (not the lane's fixed-step central difference), the
## covariance is vcov(fit, full = TRUE) selected by row name (not
## hyp_par_cov()), the map is varcorr_matrices() per block (the base
## commit's VarCorr() body), and a lane entry is matched to its
## independent counterpart by ESTIMATE VALUE, not by name.
##   Rscript dev/brmsnames-rev-varcorr-se.R      data seed 77
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
stopifnot(requireNamespace("numDeriv"))
set.seed(77)
G <- 30; m <- 12; n <- G * m
d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                g2 = factor(rep(seq_len(15), length.out = n)),
                x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                tim = factor(rep(1:4, length.out = n)))
S4 <- matrix(0.6, 4, 4); diag(S4) <- 1
Lc <- t(chol(S4 * 0.5))
U <- t(Lc %*% matrix(rnorm(4 * G), 4))
d$y <- 1 + d$x1 + U[d$g, 1] + U[d$g, 2] * d$x1 + U[d$g, 3] * d$x2 +
  U[d$g, 4] * d$x3 + rnorm(15, 0, 0.3)[d$g2] + rnorm(n)
# near +1 correlation, and a near-zero sd
u1 <- rnorm(G)
d$ynear <- 1 + u1[d$g] + 0.999 * u1[d$g] * d$x1 +
  0.001 * rnorm(G)[d$g] * d$x2 + rnorm(n)
d$y1 <- d$y; d$y2 <- 0.5 * d$y + U[d$g, 2] + rnorm(n)
d$ys <- rnorm(n, U[d$g, 1], exp(0.3 * U[d$g, 2]))
d$ytim <- rnorm(n, U[d$g, 1] + as.integer(d$tim) * 0.1)

cases <- list(
  us4_ml = list(bf(y ~ x1 + (1 + x1 + x2 + x3 | g) + (1 | g2)), FALSE),
  us4_reml = list(bf(y ~ x1 + (1 + x1 + x2 + x3 | g) + (1 | g2)), TRUE),
  near = list(bf(ynear ~ x1 + (1 + x1 + x2 | g)), FALSE),
  diag = list(bf(y ~ x1 + diag(1 + x1 + x2 | g)), FALSE),
  cs = list(bf(ytim ~ 1 + cs(tim + 0 | g)), FALSE),
  ar1 = list(bf(ytim ~ 1 + ar1(tim + 0 | g)), FALSE),
  mv_id = list(mvbf(bf(y1 ~ x1 + (1 | p | g)), bf(y2 ~ x1 + (1 | p | g)),
                    rescor = FALSE), FALSE),
  sigma_id = list(bf(ys ~ 1 + (1 | q | g), sigma ~ 1 + (1 | q | g)), FALSE)
)

ind_table <- function(fit) {
  V <- vcov(fit, full = TRUE)
  th <- fit$estimates[["theta"]]
  thn <- grep("^theta", rownames(V), value = TRUE)
  stopifnot(length(thn) == length(th))
  Vt <- V[thn, thn, drop = FALSE]
  out <- NULL
  mats0 <- varcorr_matrices(fit)
  for (b in seq_along(mats0)) {
    f <- function(t) {
      M <- varcorr_matrices(fit, t)[[b]]
      s <- sqrt(diag(M))
      C <- if (nrow(M) > 1) stats::cov2cor(M)[lower.tri(M)] else NULL
      c(s, C)
    }
    J <- numDeriv::jacobian(f, th)
    est <- f(th)
    se <- sqrt(pmax(0, rowSums((J %*% Vt) * J)))
    out <- rbind(out, data.frame(block = names(mats0)[b], est = est, se = se))
  }
  # residual sd for a scalar log-link sigma
  bn <- grep("sigma", rownames(V), value = TRUE)
  if (length(bn) == 1L && length(fit$spec$responses) == 1L) {
    e <- exp(fit$estimates$betad[1])
    out <- rbind(out, data.frame(block = "residual", est = e,
                                 se = e * sqrt(V[bn, bn])))
  }
  out
}

lane_table <- function(vc) {
  out <- NULL
  for (k in names(vc)) {
    s <- vc[[k]]$sd
    out <- rbind(out, data.frame(key = k, name = rownames(s),
                                 est = s[, "Estimate"], se = s[, "Est.Error"]))
    if (!is.null(vc[[k]]$cor)) {
      C <- vc[[k]]$cor
      K <- dim(C)[1]
      for (i in 2:K) for (j in 1:(i - 1)) {
        out <- rbind(out, data.frame(key = k,
                                     name = paste(dimnames(C)[[1]][j], dimnames(C)[[1]][i]),
                                     est = C[i, "Estimate", j], se = C[i, "Est.Error", j]))
      }
    }
  }
  out
}

for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- tryCatch(q(frm(cs[[1]], family = gaussian(), data = d, REML = cs[[2]])),
                  error = function(e) e)
  if (inherits(fit, "error")) { cat(nm, "frm ERROR", conditionMessage(fit), "\n"); next }
  vc <- tryCatch(VarCorr(fit), error = function(e) e)
  if (inherits(vc, "error")) { cat(nm, "VarCorr ERROR", conditionMessage(vc), "\n"); next }
  it <- tryCatch(ind_table(fit), error = function(e) e)
  if (inherits(it, "error")) { cat(nm, "independent ERROR", conditionMessage(it), "\n"); next }
  lt <- lane_table(vc)
  worst <- 0; unmatched <- 0; cross <- 0; rows <- 0
  for (i in seq_len(nrow(lt))) {
    if (lt$est[i] == 0 && lt$se[i] == 0) { cross <- cross + 1; next }
    hit <- which(abs(it$est - lt$est[i]) <= 1e-10 * max(1, abs(lt$est[i])))
    if (!length(hit)) { unmatched <- unmatched + 1
      cat("   unmatched:", lt$key[i], lt$name[i], lt$est[i], "\n"); next }
    rows <- rows + 1
    rel <- abs(lt$se[i] - it$se[hit[1]]) / max(it$se[hit[1]], 1e-300)
    if (rel > worst) worst <- rel
    if (rel > 1e-3) cat("   SE differs:", lt$key[i], lt$name[i], "est", lt$est[i],
                        "lane", lt$se[i], "indep", it$se[hit[1]], "rel", rel, "\n")
  }
  cat(sprintf("%-9s REML=%-5s rows %2d  cross-block zeros %2d  unmatched %d  worst rel SE diff %.3g\n",
              nm, cs[[2]], rows, cross, unmatched, worst))
}
cat("DONE\n")
