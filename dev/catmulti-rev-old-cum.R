ord_cumulative_logpmf <- function(y, eta, tau, lk, disc = 1, cs = NULL,
                                  gap_floor = FALSE) {
  if (!is.null(cs)) return(ord_cumulative_logpmf_cs(y, eta, tau, cs, lk,
                                                    disc, gap_floor))
  Fcdf <- lk$linkinv
  # Any link carrying `logit_eta` has an exact log-space difference
  # (see below), because that field turns its CDF into a logistic one.
  # cauchit is the ordinal link that does not, and does not need one:
  # its tails are polynomial, so the plain difference never saturates.
  q <- lk[["logit_eta"]]
  K1 <- length(tau)
  K <- K1 + 1L
  iK <- as.numeric(y == K)
  i1 <- as.numeric(y == 1)
  if (is.null(q)) {
    up <- Fcdf(disc * (tau[pmin(y, K1)] - eta)) * (1 - iK) + iK
    lo <- Fcdf(disc * (tau[pmax(y - 1, 1)] - eta)) * (1 - i1)
    return(log(ord_gap_floor(up - lo, gap_floor)))
  }
  # log(F(a) - F(b)) entirely in log space. `q` is the log odds of F,
  # so F is the logistic of it exactly, and for a logistic
  # F(a) - F(b) = (e^-qb - e^-qa) / ((1 + e^-qa)(1 + e^-qb)). The
  # difference of two saturated CDFs loses every digit from |eta| = 20
  # and is exactly 0 by 40; this form has no such point. For the logit
  # q is the identity and the two arguments of logspace_sub differ by
  # disc (tau_k - tau_{k-1}), which does not move with eta; for the
  # others the gap does move with eta, so the conditioning is merely
  # good rather than fixed.
  out <- i1 * log_inv_logit(q(disc * (tau[1] - eta))) +
    iK * log1m_inv_logit(q(disc * (tau[K1] - eta)))
  if (K1 >= 2L) {
    # data-only clamp into the interior categories, so the masked rows
    # still evaluate a legal (strictly ordered) threshold pair
    ym <- pmin(pmax(y, 2L), K1)
    a <- q(disc * (tau[ym] - eta))
    b <- q(disc * (tau[ym - 1L] - eta))
    out <- out + (1 - i1 - iK) * ord_log_interior(a, b, gap_floor)
  }
  out
}
