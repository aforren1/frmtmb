# Shared simulation and reference code for the HMM feasibility probe.
#
# Nothing here touches frmtmb: the numeric forward algorithm is the
# independent reference every probe is measured against.

## ---- simulation ------------------------------------------------------

#' Simulate one gaussian HMM sequence.
sim_hmm_seq <- function(Tlen, Gamma, delta, mu, sigma) {
  K <- nrow(Gamma)
  s <- integer(Tlen)
  s[1] <- sample.int(K, 1, prob = delta)
  for (t in seq_len(Tlen - 1L)) {
    s[t + 1L] <- sample.int(K, 1, prob = Gamma[s[t], ])
  }
  list(y = stats::rnorm(Tlen, mu[s], sigma[s]), state = s)
}

#' Stationary distribution of a transition probability matrix.
stat_dist <- function(Gamma) {
  K <- nrow(Gamma)
  d <- solve(t(diag(K) - Gamma) + 1, rep(1, K))
  as.vector(d)
}

## ---- numeric reference forward algorithm -----------------------------

#' Scaled forward algorithm, plain numeric. `lpmat` is T x K of
#' log emission densities.
fwd_num <- function(lpmat, Gamma, delta) {
  Tlen <- nrow(lpmat)
  a <- delta * exp(lpmat[1, ])
  s <- sum(a)
  llk <- log(s)
  a <- a / s
  for (t in seq_len(Tlen - 1L) + 1L) {
    a <- as.vector(a %*% Gamma) * exp(lpmat[t, ])
    s <- sum(a)
    llk <- llk + log(s)
    a <- a / s
  }
  llk
}

#' Log-space forward algorithm, plain numeric (independent of fwd_num:
#' different arithmetic path, so agreement is a real check).
fwd_num_log <- function(lpmat, Gamma, delta) {
  Tlen <- nrow(lpmat)
  lG <- log(Gamma)
  la <- log(delta) + lpmat[1, ]
  for (t in seq_len(Tlen - 1L) + 1L) {
    la <- apply(lG + la, 2, function(v) {
      m <- max(v)
      m + log(sum(exp(v - m)))
    }) + lpmat[t, ]
  }
  m <- max(la)
  m + log(sum(exp(la - m)))
}

#' Log emission matrix for gaussian state-dependent distributions.
lpmat_gauss <- function(y, mu, sigma) {
  vapply(seq_along(mu), function(k) stats::dnorm(y, mu[k], sigma[k],
                                                 log = TRUE),
         numeric(length(y)))
}

## ---- parameterization ------------------------------------------------

# A K x K transition matrix is K rows of a multinomial logit with the
# DIAGONAL as the reference cell (the hmmTMB/moveHMM convention): row i
# holds K - 1 free logits for the off-diagonal cells.

#' Build Gamma from a K x (K - 1) matrix of logits (numeric).
tpm_from_logits <- function(lg, K) {
  G <- matrix(0, K, K)
  for (i in seq_len(K)) {
    v <- numeric(K)
    v[-i] <- lg[i, ]
    e <- exp(v)
    G[i, ] <- e / sum(e)
  }
  G
}

#' Softmax with the first cell as reference (numeric).
softmax0 <- function(lg) {
  e <- exp(c(0, lg))
  e / sum(e)
}

## ---- RTMB forward algorithms -----------------------------------------
#
# Both take `lp` as a LIST of K advector columns (length T each), which
# sidesteps advector-matrix row extraction entirely, and `lGrows` /
# `Grows` as a list of K advector rows of the transition matrix.

#' Log-space forward recursion on the tape. `logspace_add` folded over
#' the source state keeps every step vectorized over the K target
#' states: K - 1 vector calls per time step, no elementwise assignment.
fwd_ad_log <- function(lp, lGrows, ldelta, Tlen, K) {
  la <- ldelta + vapply_ad(lp, 1L, K)
  for (t in seq_len(Tlen - 1L) + 1L) {
    acc <- lGrows[[1]] + la[1]
    if (K > 1L) {
      for (i in seq_len(K - 1L) + 1L) {
        acc <- RTMB::logspace_add(acc, lGrows[[i]] + la[i])
      }
    }
    la <- acc + vapply_ad(lp, t, K)
  }
  acc <- la[1]
  if (K > 1L) {
    for (i in seq_len(K - 1L) + 1L) acc <- RTMB::logspace_add(acc, la[i])
  }
  acc
}

#' Scaled forward recursion on the tape: the classic Zucchini scaling,
#' one division and one log per step.
fwd_ad_scale <- function(lp, Grows, delta, Tlen, K) {
  "c" <- RTMB::ADoverload("c")
  a <- delta * exp(vapply_ad(lp, 1L, K))
  s <- sum(a)
  llk <- log(s)
  a <- a / s
  for (t in seq_len(Tlen - 1L) + 1L) {
    acc <- Grows[[1]] * a[1]
    if (K > 1L) {
      for (i in seq_len(K - 1L) + 1L) acc <- acc + Grows[[i]] * a[i]
    }
    a <- acc * exp(vapply_ad(lp, t, K))
    s <- sum(a)
    llk <- llk + log(s)
    a <- a / s
  }
  llk
}

#' Row t of the emission "matrix" held as a list of K columns. `c()` is
#' the ADoverload'd one when the caller is on the tape.
vapply_ad <- function(lp, t, K) {
  "c" <- RTMB::ADoverload("c")
  out <- lp[[1]][t]
  if (K > 1L) {
    for (k in seq_len(K - 1L) + 1L) out <- c(out, lp[[k]][t])
  }
  out
}

#' Transition rows on the tape from a K x (K - 1) logit block held as a
#' plain vector in row-major order. Returns probabilities and logs.
tpm_rows_ad <- function(lgv, K) {
  "c" <- RTMB::ADoverload("c")
  Grows <- vector("list", K)
  lGrows <- vector("list", K)
  for (i in seq_len(K)) {
    idx <- (i - 1L) * (K - 1L) + seq_len(K - 1L)
    # zero in the diagonal slot, free logits elsewhere
    v <- NULL
    j2 <- 0L
    for (j in seq_len(K)) {
      el <- if (j == i) 0 else {
        j2 <- j2 + 1L
        lgv[idx[j2]]
      }
      v <- if (is.null(v)) el else c(v, el)
    }
    e <- exp(v)
    tot <- sum(e)
    Grows[[i]] <- e / tot
    lGrows[[i]] <- v - log(tot)
  }
  list(G = Grows, lG = lGrows)
}

#' Softmax over a length-(K-1) logit vector with a leading zero.
softmax0_ad <- function(lgv, K) {
  "c" <- RTMB::ADoverload("c")
  v <- 0
  if (K > 1L) for (k in seq_len(K - 1L)) v <- c(v, lgv[k])
  e <- exp(v)
  e / sum(e)
}

## ---- multiple sequences, covariate-dependent transitions -------------
#
# Convention (depmixS4's): each ROW of the transition matrix is its own
# multinomial logit with the FIRST state as the reference cell, and the
# covariate value at time t drives the transition from t to t + 1.

#' Row-wise multinomial-logit transition probabilities from a list of
#' K x K linear predictors (each a length-n numeric or advector; the
#' j = 1 entries are ignored and treated as the zero reference).
#' Returns `lg[[i]][[j]]`, log P(i -> j) at every row.
tpm_logs_ad <- function(eta, K) {
  lg <- vector("list", K)
  for (i in seq_len(K)) {
    tot <- 1  # exp(0) for the reference cell j = 1
    for (j in seq_len(K - 1L) + 1L) tot <- tot + exp(eta[[i]][[j]])
    ltot <- log(tot)
    row <- vector("list", K)
    row[[1]] <- -ltot
    for (j in seq_len(K - 1L) + 1L) row[[j]] <- eta[[i]][[j]] - ltot
    lg[[i]] <- row
  }
  lg
}

#' Log-space forward recursion over ONE sequence with time-varying
#' transitions. `rows` indexes the sequence's rows in data order; `lp`
#' and `lg` are the full-length per-state / per-cell log vectors.
fwd_ad_log_tv <- function(lp, lg, rows, ldelta, K) {
  "c" <- RTMB::ADoverload("c")
  Tl <- length(rows)
  la <- ldelta + vapply_ad(lp, rows[1], K)
  if (Tl > 1L) {
    for (s in seq_len(Tl - 1L) + 1L) {
      # covariate at the PREVIOUS row drives this step
      r <- rows[s - 1L]
      acc <- NULL
      for (i in seq_len(K)) {
        lrow <- lg[[i]][[1]][r]
        if (K > 1L) {
          for (j in seq_len(K - 1L) + 1L) lrow <- c(lrow, lg[[i]][[j]][r])
        }
        term <- lrow + la[i]
        acc <- if (is.null(acc)) term else RTMB::logspace_add(acc, term)
      }
      la <- acc + vapply_ad(lp, rows[s], K)
    }
  }
  acc <- la[1]
  if (K > 1L) {
    for (i in seq_len(K - 1L) + 1L) acc <- RTMB::logspace_add(acc, la[i])
  }
  acc
}

#' Plain-numeric multi-sequence reference with time-varying transitions.
#' `Gof(r)` returns the K x K transition matrix in force for the step
#' out of row r.
fwd_num_tv <- function(lpmat, Gof, rows, delta) {
  Tl <- length(rows)
  a <- delta * exp(lpmat[rows[1], ])
  s <- sum(a)
  llk <- log(s)
  a <- a / s
  if (Tl > 1L) {
    for (k in seq_len(Tl - 1L) + 1L) {
      a <- as.vector(a %*% Gof(rows[k - 1L])) * exp(lpmat[rows[k], ])
      s <- sum(a)
      llk <- llk + log(s)
      a <- a / s
    }
  }
  llk
}

#' Row indices of each sequence, in time order. Structure, not
#' parameters: resolved once at tape time.
hmm_seq_index <- function(gidx, tord) {
  rows <- split(seq_along(gidx), gidx)
  lapply(rows, function(r) r[order(tord[r])])
}

## ---- adaptive Gauss-Hermite quadrature -------------------------------

#' Gauss-Hermite nodes and weights for the STANDARD NORMAL weight
#' function (probabilists' convention: the weights sum to 1), by
#' Golub-Welsch on the Jacobi matrix.
gh_nodes <- function(nq) {
  i <- seq_len(nq - 1L)
  J <- matrix(0, nq, nq)
  J[cbind(i, i + 1L)] <- sqrt(i)
  J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = rev(e$values), w = rev(e$vectors[1, ]^2))
}

#' One-dimensional adaptive Gauss-Hermite quadrature of
#' `integral exp(h(b)) db`, re-centred on the mode of `h` and rescaled
#' by its curvature there. Plain GH weighted by the PRIOR fails badly
#' when the group likelihood is much sharper than the prior, which is
#' every HMM sequence longer than a handful of points.
aghq1 <- function(h, q, lo, hi) {
  op <- stats::optimize(h, c(lo, hi), maximum = TRUE, tol = 1e-10)
  bh <- op$maximum
  eps <- max(1e-4, 1e-4 * abs(bh))
  h2 <- (h(bh + eps) - 2 * op$objective + h(bh - eps)) / eps^2
  s <- if (h2 < 0) sqrt(-1 / h2) else 1
  z <- bh + s * q$x
  lv <- vapply(z, h, numeric(1)) + log(q$w) + q$x^2 / 2
  m <- max(lv)
  m + log(sum(exp(lv - m))) + log(s) + 0.5 * log(2 * pi)
}

#' Row-wise multinomial-logit transition matrices, numeric, base-1
#' reference. `B` is a K x K matrix of intercepts and `Bx` of slopes;
#' entries in column 1 are ignored.
tpm_tv_num <- function(B, Bx, x, K) {
  G <- matrix(0, K, K)
  for (i in seq_len(K)) {
    e <- c(1, exp(B[i, -1] + Bx[i, -1] * x))
    G[i, ] <- e / sum(e)
  }
  G
}
