# Short check of lane gpselimit's fix round 1: can the outward probe
# (gp_white_ridge(), the rest held at the fit) mislead?
#  P  true white fits beside another free parameter: sigma ~ z, sigma ~ x,
#     a second gp() term. Each: the gp condition, the certificate and its
#     "short" attribute, the ridge outward; any fit called short is refit
#     from lscale 1 and 3 spacings to see whether a higher maximum exists.
#  Q  the reverse: plateau starts (3e-5 of the spacing) on resolved rows
#     6 and 12 with the same extra parameters: does the held rest hide the
#     gain, so a short stop is certified?
#  C  a non-nl white fit at a doubted stop (iteration cap from a moved
#     start): the condition must be a warning, and its place in the order.
#  O  gp(by = f) row 1: the optimizer's message at code 1.
# Usage: Rscript gpselimit-short-probe.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/reviews/gpselimit-review-lib.R")
arm <- commandArgs(TRUE)[1]
rv_libs(arm)
ns <- asNamespace("frmtmb")
ts <- c(-8, -6, -4, -2, 2)
gpk <- function(conds) {
  g <- Filter(function(c) "frmtmb_gp_basis" %in% c$class, conds)
  if (!length(g)) return("none")
  paste(vapply(g, function(c) {
    paste0(substr(c$type, 1, 1), ":",
           if (grepl("stopped short", c$text)) "short" else
             if (grepl("^Boundary fit", c$text)) "boundary-msg" else "plain")
  }, ""), collapse = ";")
}
show <- function(lab, q) {
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("P", lab, arm, "ERROR", conditionMessage(fit), "\n")
    return(invisible(NULL))
  }
  ss <- rv_ses(fit)
  cat(sprintf("P %-22s %s | ll %.6f conv %d (%s) | gp [%s] | lost [%s] | conds [%s]\n",
              lab, arm, as.numeric(logLik(fit)), fit$opt$convergence,
              fit$opt$message, gpk(q$conds),
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              rv_cstr(q$conds)))
  for (bk in rv_hsgp(fit)) {
    th <- fit$estimates$theta[bk$theta_idx]
    w <- ns$gp_basis_white(bk, th)
    cert <- if (exists("gp_white_ridge", ns)) {
      tryCatch(ns$gp_white_ridge(fit, bk), error = function(e) NA)
    } else NA
    rg <- if (w) rv_ridge2(fit, bk, ts)[[1]] else rep(NA, length(ts))
    cat(sprintf("   block %-10s ls/sp %-9s white %d cert %s short %s | ridge(t=%s) %s\n",
                bk$term_label %||% "?",
                format(exp(th[2]) / bk$gp_spacing[1], digits = 3), w,
                isTRUE(cert), isTRUE(attr(cert, "short")),
                paste(ts, collapse = ","),
                paste(format(rg, digits = 3), collapse = ",")))
  }
  invisible(fit)
}
refit_far <- function(lab, fit, mk, theta_pos) {
  # from larger length scales on the white block, the rest at the default
  sp <- rv_hsgp(fit)[[1]]$gp_spacing[1]
  for (fr in c(1, 3)) {
    th <- fit$estimates$theta
    th[theta_pos] <- c(log(1), log(fr * sp))
    q <- mk(list(theta = th))
    if (inherits(q$fit, "error")) next
    cat(sprintf("   refit %s from %g sp | ll %.6f (fit %.6f, diff %.3g) ls/sp %.3g conds [%s]\n",
                lab, fr, as.numeric(logLik(q$fit)), as.numeric(logLik(fit)),
                as.numeric(logLik(q$fit)) - as.numeric(logLik(fit)),
                exp(q$fit$estimates$theta[theta_pos[2]]) / sp,
                rv_cstr(q$conds)))
  }
}
variants <- function(s, seed) {
  set.seed(seed)
  n <- nrow(s$d)
  z <- stats::runif(n, 0, 10)
  e <- stats::rnorm(n)
  list(
    sigz = list(d = transform(s$d, z = z, y = y + exp(-1.5 + 0.2 * z) * e),
                f = bf(y ~ gp(x, k = s$k), sigma ~ z)),
    sigx = list(d = transform(s$d, z = z, y = y + exp(-1.5 + 0.2 * x) * e),
                f = bf(y ~ gp(x, k = s$k), sigma ~ x)),
    gp2 = list(d = transform(s$d, z = z, y = y + sin(0.8 * z)),
               f = bf(y ~ gp(x, k = s$k) + gp(z, k = 12))))
}
for (row in c(1, 3, 21)) {
  s <- rv_sweep_row(row)
  vs <- variants(s, 7700 + row)
  for (v in names(vs)) {
    mk <- function(start = NULL) {
      rv_fit(frm(vs[[v]]$f, data = vs[[v]]$d, family = gaussian(),
                 start = start))
    }
    q <- mk()
    fit <- show(sprintf("P%02d-%s", row, v), q)
    if (is.null(fit)) next
    bk1 <- rv_hsgp(fit)[[1]]
    if (ns$gp_basis_white(bk1, fit$estimates$theta[bk1$theta_idx])) {
      refit_far(sprintf("P%02d-%s", row, v), fit, mk, bk1$theta_idx)
    }
  }
}
for (row in c(6, 12)) {
  s <- rv_sweep_row(row)
  vs <- variants(s, 7700 + row)
  for (v in names(vs)) {
    f0 <- rv_q(frm(vs[[v]]$f, data = vs[[v]]$d, family = gaussian()))
    show(sprintf("Q%02d-%s-default", row, v), list(fit = f0, conds = list()))
    bk1 <- rv_hsgp(f0)[[1]]
    th <- f0$estimates$theta
    th[bk1$theta_idx] <- c(log(30), log(3e-5 * bk1$gp_spacing[1]))
    show(sprintf("Q%02d-%s-3e-5", row, v),
         rv_fit(frm(vs[[v]]$f, data = vs[[v]]$d, family = gaussian(),
                    start = list(theta = th))))
  }
}
for (row in c(3, 21)) {
  s <- rv_sweep_row(row)
  f0 <- rv_q(frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian()))
  st <- f0$estimates[c("beta", "betad", "theta")]
  st$beta <- st$beta + 0.5
  q <- rv_fit(frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian(),
                  start = st,
                  control = frmtmb_control(optCtrl = list(iter.max = 2))))
  show(sprintf("C%02d-cap2", row), q)
  cat("   order", paste(vapply(q$conds, function(c) {
    paste0(substr(c$type, 1, 1), ":", c$class[1])
  }, ""), collapse = " > "), "\n")
}
s <- rv_sweep_row(1)
s$d$f <- factor(rep(c("a", "b"), length.out = nrow(s$d)))
show("O01-by", rv_fit(frm(y ~ gp(x, by = f, k = s$k), data = s$d)))
