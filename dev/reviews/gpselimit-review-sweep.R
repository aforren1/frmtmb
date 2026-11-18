# Review of lane gpselimit, checks 1, 3 and 4: the gpfield review's 30
# gaussian approximate gp() fits (dev/reviews/gpfield-review-basis.R, same
# data), read through the API (sdr_of(), vcov(), summary()) and not from
# warning text. Per row: every condition, the SE of every outer parameter,
# the lost reasons and kinds, the review's own "se_lost" metric (warning
# text), the objective along the ridge, and on white rows the product
# 2 log sdgp + log lscale through hypothesis() against a numerical Hessian
# in (psi, ridge) coordinates.
# Usage: Rscript gpselimit-review-sweep.R <lane|base> [jitter]
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/reviews/gpselimit-review-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]
jit <- if (length(a) >= 2) as.numeric(a[2]) else 0
rv_libs(arm)
ns <- asNamespace("frmtmb")
ts <- c(-4, -2, -1, 1, 2, 4, 8, 12, 20)
res <- list()
for (row in 1:30) {
  s <- rv_sweep_row(row, jit)
  q <- rv_fit(frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian()))
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("RVSWEEP", row, "ERROR", conditionMessage(fit), "\n"); next
  }
  bk <- rv_hsgp(fit)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  sp <- bk$gp_spacing
  ss <- rv_ses(fit)
  v <- rv_q(vcov(fit))
  sm <- rv_q(capture.output(print(summary(fit))))
  foot <- grep("white limit|boundary|singular", sm, value = TRUE)
  kind <- fit$cache$se_boundary_kind
  review_metric <- any(grepl("Standard errors are not available",
                             vapply(q$conds, `[[`, "", "text")))
  white <- ns$gp_basis_white(bk, fit$estimates$theta[bk$theta_idx])
  rg <- rv_ridge2(fit, bk, ts)[[1]]
  hyp <- c(NA, NA); hess_se <- NA
  if (white) {
    vn <- ns$gp_brms_terms(fit)[[1]]
    h <- sprintf("2 * log(%s) + log(%s) = 0", vn$sdgp[1], vn$lscale[1, 1])
    hr <- tryCatch(rv_q(hypothesis(fit, h, class = ""))$hypothesis,
                   error = function(e) NULL)
    if (!is.null(hr)) hyp <- c(hr$Estimate, hr$Est.Error)
    om <- ns$outer_par_map(fit)
    pos <- which(om$comp == "theta")[bk$theta_idx]
    p <- fit$opt$par
    H <- stats::optimHess(p, fit$obj$fn, fit$obj$gr)
    A <- diag(length(p))
    A[, pos[1]] <- 0; A[pos[1], pos[1]] <- 0.5
    A[, pos[2]] <- 0; A[pos[1], pos[2]] <- 0.5; A[pos[2], pos[2]] <- -1
    Hn <- t(A) %*% H %*% A
    keep <- setdiff(seq_along(p), pos[2])
    Vn <- solve(Hn[keep, keep])
    hess_se <- sqrt(Vn[match(pos[1], keep), match(pos[1], keep)])
    fit$obj$fn(p)
  }
  res[[row]] <- list(row = row, est = fit$estimates, par = fit$opt$par,
                     ll = as.numeric(logLik(fit)), se = ss$se,
                     lost = ss$lost, kind = kind, vcov = v,
                     conds = q$conds, foot = foot, ridge = rg, hyp = hyp,
                     hess_se = hess_se)
  cat(sprintf(paste0("RVSWEEP %2d %s j%g k %2d | ll %.6f conv %d | ",
                     "ls/sp %.3g sd %.4g white %d | se %s | lost [%s] | ",
                     "kind [%s] | conds [%s] | review_metric %d | ",
                     "ridge(t=%s) %s | hyp %s se %s hessSE %s | foot %s\n"),
              row, arm, jit, s$k, as.numeric(logLik(fit)),
              fit$opt$convergence, exp(th[2]) / sp, exp(th[1]), white,
              paste(format(ss$se, digits = 4), collapse = ","),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              paste(paste0(names(kind), "=", kind), collapse = ","),
              rv_cstr(q$conds), review_metric,
              paste(ts, collapse = ","),
              paste(format(rg, digits = 3), collapse = ","),
              format(hyp[1], digits = 6), format(hyp[2], digits = 6),
              format(hess_se, digits = 6),
              paste(substr(foot, 1, 70), collapse = " / ")))
}
saveRDS(res, file.path(RV_OUT, sprintf("sweep-%s-j%g.rds", arm, jit)))
