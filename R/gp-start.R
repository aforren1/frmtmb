#' Where the descent of a fit with an exact gp() begins.
#'
#' A gaussian response integrates the field out in closed form, so the marginal
#' of the model is cheap on numbers and on a small tape. The fit reads it to
#' choose a start and nothing more: `fit$obj` stays the z tape, every reader and
#' the maximum stay the tape's. From `d / 4` or a start read at unit sd a
#' descent can end at a long field that took the wiggle into sigma, 417 nats
#' short with conv 0 on `2 sin(3x) + x`, n 80, noise 1e-3
#' (dev/ncgp-findings.md). Everything the marginal does not cover gets
#' gp_start_scan_units().
#'
#' @name gp-start
#' @noRd
NULL

# Above this many distinct rows of the random design, the marginal's dense
# solve costs more than the descent it would save (dev/ncgp-findings.md, cost).
gp_marginal_max_rows <- 1000L

gp_marginal_blocks <- c("gp", "us", "diag", "cs", "homdiag", "smooth", "hsgp")

#' The one predictor of mu, with sigma a constant, when the marginal of the
#' model is a multivariate normal; otherwise `NULL`.
#'
#' @noRd
gp_marginal_gate <- function(frame) {
  spec <- frame[["spec"]]
  if (length(spec$responses) != 1L) return(NULL)
  resp <- spec$responses[[1L]]
  fam <- resp$family
  if (!identical(fam[["family"]], "gaussian")) return(NULL)
  if (length(resp$nlpars) || length(frame[["autocor"]]) ||
        length(frame[["me"]]$terms) || length(frame[["extra_names"]])) {
    return(NULL)
  }
  fixed <- length(frame[["betad_fixed_idx"]]) > 0L
  if (length(setdiff(names(frame[["map"]]), if (fixed) "betad"))) return(NULL)
  tpl <- frame[["par_template"]]
  if (!all(names(tpl) %in% c("beta", "betad", "b", "theta"))) return(NULL)
  aterms <- frame[["aterm_values"]][[1L]]
  if (!all(names(aterms) %in% c("se", "se_sigma", "weights"))) return(NULL)
  ok <- vapply(frame[["re_blocks"]], function(bk) {
    bk[["covstruct"]] %in% gp_marginal_blocks
  }, NA)
  if (!all(ok)) return(NULL)
  lps <- frame[["linpreds"]]
  mu <- Filter(function(l) identical(l[["dpar"]], "mu"), lps)
  sg <- Filter(function(l) identical(l[["dpar"]], "sigma"), lps)
  if (length(mu) != 1L || length(lps) != length(mu) + length(sg)) return(NULL)
  mu <- mu[[1L]]
  if (!is.null(mu[["constant"]]) || !identical(mu[["par"]], "beta")) {
    return(NULL)
  }
  if (!isTRUE(all.equal(mu[["link"]]$linkfun(c(2, 3)), c(2, 3)))) return(NULL)
  p <- length(mu[["idx"]])
  if (p != NCOL(mu[["X"]]) || p != length(tpl[["beta"]])) return(NULL)
  has_sigma <- FALSE
  if (length(sg) && !fixed) {
    sg <- sg[[1L]]
    if (NCOL(sg[["X"]]) != 1L || length(tpl[["betad"]]) != 1L ||
          !is.null(sg[["constant"]]) || any(sg[["offset"]] != 0)) {
      return(NULL)
    }
    has_sigma <- TRUE
  }
  list(mu = mu, has_sigma = has_sigma, aterms = aterms)
}

#' The model's marginal over the distinct rows of the random design.
#'
#' Rows with one row of `Z`, one `se` and one weight are sufficient through
#' their mean and a within sum of squares, so the matrix to factor has the size
#' of the positions and not of the rows: 1000 rows on 500 positions factor 500.
#'
#' `V(theta, noise)` is the marginal covariance of the means, on numbers or on
#' the tape; `noise` is sigma^2 plus the row's `se^2` handled inside.
#'
#' @noRd
gp_marginal_setup <- function(frame, max_rows = gp_marginal_max_rows) {
  g <- gp_marginal_gate(frame)
  if (is.null(g)) return(NULL)
  mu <- g$mu
  y <- as.numeric(frame[["y"]][[1L]]) - as.numeric(mu[["offset"]] %||% 0)
  X <- as.matrix(mu[["X"]])
  Z <- mu[["Z"]]
  n <- length(y)
  if (!all(is.finite(y)) || NROW(X) != n || NROW(Z) != n) return(NULL)
  se2 <- g$aterms[["se"]]
  se2 <- if (is.null(se2)) rep(0, n) else as.numeric(se2)^2
  w <- g$aterms[["weights"]]
  w <- if (is.null(w)) rep(1, n) else as.numeric(w)
  gp_marginal_data(y, X, Z, frame[["re_blocks"]], se2, w, g$has_sigma,
                   mu[["idx"]], frame, max_rows)
}

#' The marginal's sufficient statistics and blocks from `y`, `X`, `Z`.
#'
#' @noRd
gp_marginal_data <- function(y, X, Z, re_blocks, se2, w, has_sigma, idx, frame,
                             max_rows = gp_marginal_max_rows) {
  n <- length(y)
  # a fixed vector, not a draw, so a fit leaves the RNG alone; two rows of Z
  # read one number only when they are the same row
  if (all(se2 == se2[1L]) && all(w == w[1L])) {
    h <- as.vector(Z %*% cos(seq_len(ncol(Z))))
    J <- match(h, unique(h))
  } else {
    J <- seq_len(n)
  }
  u <- max(J)
  if (u > max_rows) return(NULL)
  first <- match(seq_len(u), J)
  nj <- tabulate(J, u)
  ybar <- rowsum(y, J, reorder = TRUE)[, 1L] / nj
  Xbar <- if (ncol(X)) rowsum(X, J, reorder = TRUE) / nj else matrix(0, u, 0L)
  yw <- y - ybar[J]
  Xw <- X - Xbar[J, , drop = FALSE]
  Zu <- Z[first, , drop = FALSE]
  blocks <- lapply(re_blocks, function(bk) {
    Zb <- Zu[, bk[["c_idx"]], drop = FALSE]
    if (length(bk[["theta_idx"]]) == 1L &&
          (bk[["dim"]] == 1L || identical(bk[["covstruct"]], "smooth"))) {
      return(list(kind = "scale", idx = bk[["theta_idx"]],
                  ZZ = as.matrix(Matrix::tcrossprod(Zb))))
    }
    if (identical(bk[["covstruct"]], "gp")) {
      # a row reads one position, times a by-variable: an index, not a product
      nz <- Matrix::summary(methods::as(Zb, "TsparseMatrix"))
      if (nrow(nz) == u && all(tabulate(nz$i, u) == 1L)) {
        o <- order(nz$i)
        return(list(kind = "gp", bk = bk, pos = nz$j[o],
                    mm = tcrossprod(nz$x[o])))
      }
    }
    list(kind = "general", bk = bk, Zb = as.matrix(Zb))
  })
  list(y = ybar, X = Xbar, n = n, u = u, nj = nj, w = w[first],
       se2 = se2[first], ywyw = sum(yw^2), Xwyw = crossprod(Xw, yw),
       XwXw = crossprod(Xw), logn = sum(log(nj)), blocks = blocks,
       has_sigma = has_sigma, idx = idx, frame = frame)
}

#' The field and random-effect part of the marginal covariance of the means.
#'
#' @noRd
gp_marginal_cov <- function(s, theta) {
  "[<-" <- RTMB::ADoverload("[<-")
  V <- 0
  for (b in s$blocks) {
    if (identical(b$kind, "scale")) {
      V <- V + exp(2 * theta[b$idx]) * b$ZZ
    } else if (identical(b$kind, "gp")) {
      C <- cs_entry("gp")$vcov(theta[b$bk[["theta_idx"]]], b$bk)
      V <- V + C[b$pos, b$pos, drop = FALSE] * b$mm
    } else {
      bk <- b$bk
      C <- cs_entry(bk[["covstruct"]])$vcov(theta[bk[["theta_idx"]]], bk)
      C <- RTMB::matrix(C, bk[["dim"]], bk[["dim"]])
      for (l in seq_len(bk[["n_levels"]])) {
        ii <- (l - 1L) * bk[["dim"]] + seq_len(bk[["dim"]])
        Zl <- b$Zb[, ii, drop = FALSE]
        V <- V + Zl %*% C %*% t(Zl)
      }
    }
  }
  V
}

#' The marginal's negative log density, on the tape. Same constant as the fit's
#' objective, so the two read one number at one point on a gaussian response.
#'
#' @noRd
gp_marginal_nll <- function(s) {
  function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    sig2 <- if (s$has_sigma) exp(2 * p$betad) else 0
    noise <- (sig2 / s$w + s$se2) / s$nj
    V <- gp_marginal_cov(s, p$theta) + RTMB::diag(noise, s$u)
    e <- if (ncol(s$X)) s$y - as.vector(s$X %*% p$beta) else s$y
    out <- -sum(RTMB::dmvnorm(e, 0, V, log = TRUE))
    if (s$n > s$u) {
      within <- s$ywyw
      if (ncol(s$X)) {
        within <- within - 2 * sum(p$beta * s$Xwyw) +
          sum(p$beta * as.vector(s$XwXw %*% p$beta))
      }
      den <- sig2 / s$w[1L] + s$se2[1L]
      out <- out + 0.5 * within / den +
        0.5 * (s$n - s$u) * log(2 * pi * den) + 0.5 * s$logn
    }
    out
  }
}

#' One cell of the coarse scan on numbers: beta by generalized least squares, and
#' with `sigma` free and no `se()` a common scale of the whole covariance read in
#' closed form, so the grid has only the length scales and the noise ratio.
#'
#' `theta` holds the sds on the scale `1` when `lr` is given. With a known noise
#' `noise` is sigma^2 and the scale `s` is the caller's. Returns the value and the
#' point.
#'
#' @noRd
gp_marginal_cell <- function(s, theta, lr, s0 = NULL) {
  r <- exp(lr)
  known <- !is.null(s0)
  Cb <- gp_marginal_cov(s, theta)
  den <- if (!known) r / s$w else
    (if (s$has_sigma) r * s0^2 else 0) / s$w + s$se2
  V <- (if (known) s0^2 * Cb else Cb)
  diag(V) <- diag(V) + den / s$nj
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(list(value = Inf))
  wy <- backsolve(R, s$y, transpose = TRUE)
  p <- ncol(s$X)
  wX <- if (p) backsolve(R, s$X, transpose = TRUE) else matrix(0, s$u, 0L)
  den1 <- den[1L]
  has_within <- s$n > s$u
  XtX <- crossprod(wX) + (if (has_within) s$XwXw / den1 else 0)
  Xty <- crossprod(wX, wy) + (if (has_within) s$Xwyw / den1 else 0)
  beta <- if (p) tryCatch(solve(XtX, Xty), error = function(e) NULL) else numeric(0)
  if (is.null(beta) || !all(is.finite(beta))) return(list(value = Inf))
  we <- if (p) wy - wX %*% beta else wy
  within <- if (has_within) {
    s$ywyw - (if (p) 2 * sum(beta * s$Xwyw) - sum(beta * (s$XwXw %*% beta)) else 0)
  } else 0
  q <- sum(we^2) + within / (if (has_within) den1 else 1)
  ldet <- 2 * sum(log(diag(R)))
  if (has_within) ldet <- ldet + (s$n - s$u) * log(den1) + s$logn
  if (known) {
    v <- 0.5 * (s$n * log(2 * pi) + ldet + q)
    return(list(value = v, beta = drop(beta), s = s0))
  }
  s2 <- q / s$n
  if (!isTRUE(s2 > 0)) return(list(value = Inf))
  v <- 0.5 * (s$n * log(2 * pi * s2) + ldet + s$n)
  list(value = v, beta = drop(beta), s = sqrt(s2))
}

#' A theta with each exact gp()'s sd at `lsd`, each other block's sds one step
#' below it, and every length scale where `rho` says.
#'
#' @noRd
gp_marginal_theta <- function(frame, lsd, rho = NULL) {
  th <- frame[["par_template"]][["theta"]]
  for (i in seq_along(frame[["re_blocks"]])) {
    bk <- frame[["re_blocks"]][[i]]
    if (identical(bk[["covstruct"]], "gp")) {
      th[bk[["theta_idx"]][1L]] <- lsd
      r <- rho[[as.character(i)]]
      if (is.null(r)) {
        r <- sqrt(max(Reduce(`+`, bk[["aux_D2"]]))) / 4
      }
      th[bk[["theta_idx"]][-1L]] <- log(r)
    } else {
      si <- cs_entry(bk[["covstruct"]])[["sd_idx"]]
      si <- if (is.function(si)) si(bk) else si
      if (length(si)) th[bk[["theta_idx"]][si]] <- lsd - 1
    }
  }
  th
}

#' The marginal's maximum from a coarse grid: local maxima from the best cell of
#' `distinct` different length scales (six for the warm start, 1 of 744 rows
#' short, and four for the working start, where six lose on low-count Poisson and
#' binomial) and from the best `keep` cells overall, the
#' best kept. The best cells overall usually share a basin, and on fresh curves
#' that alone stopped up to 40 nats short (dev/ncgp-findings.md, round 2). The grid
#' sweeps one block's length scale at a time with the others at the best so far.
#' A list `value`, `beta`, `betad`, `theta`, or `NULL` where the marginal does not
#' apply or a step fails.
#'
#' @noRd
gp_marginal_max <- function(frame, keep = 3L, distinct = 6L) {
  s <- tryCatch(gp_marginal_setup(frame), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  tryCatch(gp_marginal_opt(s, keep, distinct), error = function(e) NULL)
}

gp_marginal_opt <- function(s, keep = 3L, distinct = 4L) {
  frame <- s$frame
  blocks <- frame[["re_blocks"]]
  gpi <- which(vapply(blocks, function(b) {
    identical(b[["covstruct"]], "gp")
  }, NA))
  has_se <- any(s$se2 > 0)
  lrs <- log(10^c(1, -1, -3, -6, -9))
  res <- if (ncol(s$X)) {
    s$y - s$X %*% qr.solve(s$X * sqrt(s$nj), s$y * sqrt(s$nj))
  } else s$y
  sy <- sqrt(sum(s$nj * res^2) / s$n)
  if (!isTRUE(sy > 0)) sy <- 1
  scales <- if (has_se) sy * c(0.3, 1, 3) else list(NULL)
  rho <- list()
  cells <- list()
  for (i in gpi) {
    d <- sqrt(max(Reduce(`+`, blocks[[i]][["aux_D2"]])))
    if (!isTRUE(d > 0)) next
    best <- list(value = Inf)
    for (r in d / 2^(0:10)) {
      rho2 <- rho
      rho2[[as.character(i)]] <- r
      for (s0 in scales) {
        for (lr in lrs) {
          lsd <- if (is.null(s0)) 0 else log(s0)
          th <- gp_marginal_theta(frame, lsd, rho2)
          cv <- gp_marginal_cell(s, th, lr, s0)
          if (!is.finite(cv$value)) next
          cells[[length(cells) + 1L]] <- list(value = cv$value, rho = rho2,
                                              lr = lr, beta = cv$beta, s = cv$s)
          if (cv$value < best$value) best <- list(value = cv$value, r = r)
        }
      }
    }
    if (is.finite(best$value)) rho[[as.character(i)]] <- best$r
  }
  if (!length(cells)) return(NULL)
  vals <- vapply(cells, `[[`, 0, "value")
  rk <- vapply(cells, function(cc) {
    paste(signif(unlist(cc$rho), 6), collapse = "_")
  }, "")
  o0 <- order(vals)
  o0 <- o0[!duplicated(rk[o0])]
  ord <- unique(c(o0[seq_len(min(distinct, length(o0)))],
                  order(vals)[seq_len(min(keep, length(vals)))]))
  nll <- gp_marginal_nll(s)
  old <- RTMB::TapeConfig()
  on.exit(RTMB::TapeConfig(old), add = TRUE)
  RTMB::TapeConfig(matmul = "atomic")
  out <- NULL
  for (k in ord) {
    cc <- cells[[k]]
    p <- list(theta = gp_marginal_theta(frame, log(cc$s), cc$rho))
    if (ncol(s$X)) p$beta <- cc$beta
    if (s$has_sigma) p$betad <- log(cc$s) + cc$lr / 2
    ob <- RTMB::MakeADFun(nll, p, silent = TRUE)
    # beta in y's units, as everything else here, so the descent does not
    # depend on them
    sc <- ifelse(names(ob$par) == "beta", 1 / sy, 1)
    ctl <- list(iter.max = 300L, eval.max = 400L)
    o <- tryCatch(stats::nlminb(ob$par, ob$fn, ob$gr, scale = sc, control = ctl),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$objective)) next
    if (is.null(out) || o$objective < out$value) {
      pl <- ob$env$parList(o$par)
      out <- list(value = o$objective, beta = as.numeric(pl$beta),
                  betad = pl$betad, theta = as.numeric(pl$theta))
    }
  }
  if (is.null(out)) return(NULL)
  if (s$has_sigma && !has_se) out <- gp_marginal_polish(s, out)
  bt <- numeric(length(frame[["par_template"]][["beta"]]))
  bt[s$idx] <- out$beta
  out$beta <- bt
  out
}

#' A last descent with `beta` and the common scale read in closed form at each
#' step. On a trend the field also explains, `exp(0.8 x)`, the intercept and the sd
#' move along a long ridge, 38784 against 371, which a descent that carries them
#' stops on: 0.67 nats short with no second point there (dev/ncgp-findings.md).
#'
#' @noRd
gp_marginal_polish <- function(s, out) {
  blocks <- s$frame[["re_blocks"]]
  sds <- unlist(lapply(blocks, function(bk) {
    si <- cs_entry(bk[["covstruct"]])[["sd_idx"]]
    bk[["theta_idx"]][if (is.function(si)) si(bk) else si]
  }))
  i0 <- blocks[[which(vapply(blocks, function(b) {
    identical(b[["covstruct"]], "gp")
  }, NA))[1L]]][["theta_idx"]][1L]
  shift <- out$theta[i0]
  th0 <- out$theta
  th0[sds] <- th0[sds] - shift
  free <- setdiff(seq_along(th0), i0)
  f <- function(q) {
    th <- th0
    th[free] <- q[seq_along(free)]
    gp_marginal_cell(s, th, q[length(q)])$value
  }
  q0 <- c(th0[free], 2 * (out$betad - shift))
  ctl <- list(iter.max = 60L, eval.max = 150L, rel.tol = 1e-12)
  o <- tryCatch(stats::nlminb(q0, f, control = ctl),
                error = function(e) NULL)
  if (is.null(o) || !is.finite(o$objective)) return(out)
  # along that ridge the value rounds at 1e-3, and a gradient by differences reads
  # the rounding: nlminb stopped 0.2 nats short where a simplex reaches the maximum
  if (o$convergence != 0L || grepl("X-conv", o$message, fixed = TRUE)) {
    nm <- tryCatch(stats::optim(o$par, f, method = "Nelder-Mead",
                                control = list(maxit = 400L, reltol = 1e-14)),
                   error = function(e) NULL)
    if (!is.null(nm) && is.finite(nm$value) && nm$value < o$objective) {
      o$par <- nm$par
      o$objective <- nm$value
    }
  }
  if (o$objective >= out$value) return(out)
  th <- th0
  th[free] <- o$par[seq_along(free)]
  cv <- gp_marginal_cell(s, th, o$par[length(o$par)])
  th[sds] <- th[sds] + log(cv$s)
  list(value = cv$value, beta = cv$beta, theta = th,
       betad = array(log(cv$s) + o$par[length(o$par)] / 2))
}

#' The outer gradient of a fit whose objective is the marginal of an exact
#' gp(), from the marginal's own tape, or `NULL`. The z tape's gradient runs an
#' inner solve over a dense Hessian each time: at 500 positions the outer
#' Hessian's differences took 11.3 s of a 25 s fit, and the marginal's tape
#' 0.87 s (dev/reviews/2026-11-03-ncgp.md, cost). The fit's own maximum
#' likelihood only, with no prior: anything else is another objective. The two
#' must read one value at `at`, or the z tape keeps the job.
#'
#' The marginal's gradient rounds with the condition of its covariance `V`,
#' where the z tape's does not, so `V` must be well conditioned too: its bound
#' `(trace of the field part + largest noise) / smallest noise` at most
#' `gp_marginal_kappa_max`.
#'
#' @noRd
gp_marginal_outer_gr <- function(fit, at) {
  fr <- fit$frame
  if (!any(vapply(fr[["re_blocks"]] %||% list(), function(bk) {
    identical(bk[["covstruct"]], "gp")
  }, NA))) {
    return(NULL)
  }
  if (isTRUE(fit$REML) || isTRUE(fit$quadrature) || !is.null(fit$importance) ||
        !is.null(fit$beta_whiten) || isTRUE(fit$control$profile) ||
        NROW(fit$prior) > 0L) {
    return(NULL)
  }
  pn <- names(at)
  if (!length(pn) || !all(pn %in% c("beta", "betad", "theta"))) return(NULL)
  s <- tryCatch(gp_marginal_setup(fr, max_rows = Inf),
                error = function(e) NULL)
  if (is.null(s) || s$has_sigma != ("betad" %in% pn)) return(NULL)
  sig2 <- if (s$has_sigma) exp(2 * at[[which(pn == "betad")]]) else 0
  noise <- (sig2 / s$w + s$se2) / s$nj
  th <- unname(at[pn == "theta"])
  vf <- tryCatch(diag(as.matrix(gp_marginal_cov(s, th))),
                 error = function(e) NA_real_)
  kap <- (sum(vf) + max(noise)) / min(noise)
  if (!isTRUE(kap <= gp_marginal_kappa_max)) return(NULL)
  nll <- gp_marginal_nll(s)
  idx <- s$idx
  f <- function(p) {
    if (!is.null(p$beta)) p$beta <- p$beta[idx]
    nll(p)
  }
  pl <- lapply(c(beta = "beta", betad = "betad", theta = "theta"),
               function(k) unname(at[pn == k]))
  pl <- pl[lengths(pl) > 0L]
  old <- RTMB::TapeConfig()
  ob <- tryCatch({
    RTMB::TapeConfig(matmul = "atomic")
    RTMB::MakeADFun(f, pl, silent = TRUE)
  }, error = function(e) NULL, finally = RTMB::TapeConfig(old))
  if (is.null(ob) || !identical(names(ob$par), pn)) return(NULL)
  v <- tryCatch(ob$fn(at), error = function(e) NA_real_)
  v0 <- tryCatch(fit$obj$fn(at), error = function(e) NA_real_)
  if (!isTRUE(abs(v - v0) <= gp_marginal_value_tol * max(1, abs(v0)))) {
    return(NULL)
  }
  function(q) ob$gr(q)
}

# How far apart, relative to the value, the marginal and the z tape may read at
# one point before the z tape keeps the outer Hessian. The inner problem is
# gaussian and linear, so they differ by the inner solve's rounding only.
gp_marginal_value_tol <- 1e-9

# The largest bound on the condition of the marginal's covariance at which it
# gives the outer Hessian. Against the z tape on 377 fits of the ncgp review's
# sets, the standard errors moved at most 8.7e-8 below 1e6, 5.6e-4 below 1e9,
# and beyond 1e9 by up to 1.0, with 5 sets of NaN standard errors changed; at
# noise 1e-3 the marginal's were 35 percent off an exact 200-bit Hessian and
# the z tape's 3e-5 (dev/fitedge-findings.md, item 3).
gp_marginal_kappa_max <- 1e6

#' The value and length scale a link-scale response reads for each exact gp() a
#' response's mean holds, for what the marginal's gate refuses: another family, a
#' non-identity link, `sigma ~`, `ar()`, `mi()`, `me()`, `mo()`, several responses.
#' Only the blocks' `theta` is taken. A start from where the objective reads unit
#' sd stopped up to 203 nats short on those, conv 0 and silent; 0.78.0 stopped short
#' on the same rows, and on `mi()`, `me()` and `mo()` it did not
#' (dev/ncgp-findings.md, round 2).
#'
#' A block in a parameter that a nonlinear body of the mean reads gets the same
#' start, from gp_body_start(): with no start for it, `a + b * z` with
#' `a ~ gp(x)` stopped 100 to 114 nats short, conv 0 and silent
#' (dev/fitedge-findings.md, item 1).
#'
#' `tpl` is the template the descent starts from, which places the body's other
#' coefficients. `body = FALSE` leaves those blocks out, as a caller's `start`
#' does: from a far start the body read a field of sd 1e-11 and stopped 26 and
#' 125 nats below 0.79.1 (dev/reviews/2026-11-06-fitedge.md, B1).
#'
#' @return `list(idx, theta, body)`, the positions in the template's `theta`,
#'   their values, and which of the positions the body start set, or `NULL`.
#' @noRd
gp_start_working <- function(frame, tpl = frame[["par_template"]],
                             body = TRUE) {
  blocks <- Filter(function(bk) identical(bk[["covstruct"]], "gp"),
                   frame[["re_blocks"]] %||% list())
  keys <- vapply(blocks, gp_block_lp, "", frame = frame)
  idx <- integer(0)
  theta <- numeric(0)
  from_body <- integer(0)
  for (key in setdiff(unique(keys), "")) {
    lp <- frame[["linpreds"]][[key]]
    if (!is.null(lp[["nl_body"]]) || !is.null(lp[["constant"]])) next
    rb <- blocks[keys == key]
    at <- 0L
    from <- integer(0)
    for (i in seq_along(rb)) {
      ti <- rb[[i]][["theta_idx"]]
      rb[[i]][["theta_idx"]] <- at + seq_along(ti)
      at <- at + length(ti)
      from <- c(from, ti)
    }
    fr <- list(re_blocks = rb, par_template = list(theta = numeric(at), beta = 0))
    if (!identical(lp[["dpar"]], "mu")) {
      if (!body) next
      m <- tryCatch(gp_body_start(frame, key, tpl, rb, fr),
                    error = function(e) NULL)
      if (is.null(m) || !all(is.finite(m$theta))) next
      idx <- c(idx, from)
      theta <- c(theta, m$theta)
      from_body <- c(from_body, from)
      next
    }
    yw <- gp_working_y(frame, lp)
    if (is.null(yw)) next
    off <- as.numeric(lp[["offset"]] %||% 0)
    yw <- yw - off
    ok <- is.finite(yw)
    if (sum(ok) < 3L || NROW(lp[["Z"]]) != length(yw)) next
    n <- sum(ok)
    m <- tryCatch(gp_marginal_opt(
      gp_marginal_data(yw[ok], matrix(1, n, 1L), lp[["Z"]][ok, , drop = FALSE],
                       rb, rep(0, n), rep(1, n), TRUE, 1L, fr)),
      error = function(e) NULL)
    if (is.null(m) || !all(is.finite(m$theta))) next
    idx <- c(idx, from)
    theta <- c(theta, m$theta)
  }
  if (length(idx)) list(idx = idx, theta = theta, body = from_body)
}

#' The marginal's maximum for the blocks `rb` of a parameter `a` that a
#' nonlinear body of the mean reads, or `NULL`: gp_marginal_opt() on the body
#' linearized at the template, one Gauss-Newton step. With `y_w` the mean's
#' working response and `eta0` the body at the template with every field at 0,
#' `y_w - eta0 + J_a a0 + sum(J_k b_k0)` is linear in `a`'s coefficients, in
#' each other coefficient `b_k` the body reads and in the field, with the
#' field's rows scaled by `J_a`; the derivatives are central differences. On
#' `a + b * z` that is the linear spelling's marginal, so the body starts where
#' `y ~ z + gp(x)` does: a start from `y` alone read 0.44 nats short on one of
#' 24 rows (dev/fitedge-findings.md, item 1). With the second descent from
#' 0.79.1's start beside it (gp_second_descent()), a second step from the
#' first step's coefficients never read higher on 48 rows of `exp(a)`,
#' `a * exp(b * z)` and the logistic body, and read up to 44 nats lower on 4;
#' one step alone, with no second descent, is worse than 0.79.1 on 6 of 96
#' rows (dev/fitedge-findings.md, fix round). Rows where the body barely reads
#' `a` (`|J_a|` under 1e-3 of its largest) are left out, and a coefficient the
#' design cannot tell from the others keeps its start.
#'
#' @noRd
gp_body_start <- function(frame, key, tpl, rb, fr) {
  lp <- frame[["linpreds"]][[key]]
  rn <- lp[["resp"]]
  nl <- frame[["spec"]]$responses[[rn]]$nlpars %||% character(0)
  body <- frame[["linpreds"]][[linpred_key(rn, "mu")]]
  if (!lp[["dpar"]] %in% nl || is.null(body) || is.null(body[["nl_body"]]) ||
        !identical(body[["dpar"]], "mu") ||
        NROW(lp[["Z"]]) != NROW(lp[["X"]])) {
    return(NULL)
  }
  ym <- gp_working_y(frame, body)
  if (is.null(ym)) return(NULL)
  off <- as.numeric(lp[["offset"]] %||% 0)
  # the body with every row of `a` shifted by h
  ev <- function(t, h = 0) {
    f2 <- frame
    f2[["linpreds"]][[key]][["offset"]] <- off + h
    e <- tryCatch(suppressWarnings(nl_body_at(f2, body, t)),
                  error = function(e) NULL)
    if (length(e) == length(ym)) as.numeric(e)
  }
  own <- lp[["idx"]]
  bc <- nl_body_coefs(frame, body)$coefs
  others <- do.call(rbind, lapply(names(bc), function(par) {
    pos <- bc[[par]]
    if (identical(par, lp[["par"]])) pos <- setdiff(pos, own)
    if (length(pos)) data.frame(par = par, pos = pos, stringsAsFactors = FALSE)
  }))
  Xa <- as.matrix(lp[["X"]])
  ba <- tpl[[lp[["par"]]]][own]
  a0 <- off + if (length(ba)) as.numeric(Xa %*% ba) else 0
  h <- 1e-4 * max(1, abs(a0))
  e0 <- ev(tpl)
  e1 <- ev(tpl, h)
  e2 <- ev(tpl, -h)
  if (is.null(e0) || is.null(e1) || is.null(e2)) return(NULL)
  Ja <- (e1 - e2) / (2 * h)
  yt <- ym - e0 + Ja * (a0 - off)
  Xt <- if (ncol(Xa)) Ja * Xa else matrix(0, length(Ja), 0L)
  for (k in seq_len(NROW(others))) {
    par <- others$par[k]
    v <- tpl[[par]][others$pos[k]]
    hk <- 1e-4 * max(1, abs(v))
    tp <- tm <- tpl
    tp[[par]][others$pos[k]] <- v + hk
    tm[[par]][others$pos[k]] <- v - hk
    ep <- ev(tp)
    em <- ev(tm)
    if (is.null(ep) || is.null(em)) next
    Jk <- (ep - em) / (2 * hk)
    if (!all(is.finite(Jk))) next
    yt <- yt + Jk * v
    Xt <- cbind(Xt, Jk)
  }
  ok <- is.finite(yt) & is.finite(Ja) & apply(is.finite(Xt), 1L, all)
  if (sum(ok) < 3L) return(NULL)
  ok <- ok & abs(Ja) > 1e-3 * max(abs(Ja[ok]))
  n <- sum(ok)
  if (n < 3L) return(NULL)
  # a column the rows cannot tell from the others, `a + b` with two
  # intercepts, keeps its start
  q <- qr(Xt[ok, , drop = FALSE])
  keep <- sort(q$pivot[seq_len(q$rank)])
  Zt <- Matrix::Diagonal(x = Ja[ok]) %*% lp[["Z"]][ok, , drop = FALSE]
  f2 <- fr
  f2[["par_template"]][["beta"]] <- numeric(length(keep))
  m <- tryCatch(gp_marginal_opt(
    gp_marginal_data(yt[ok], Xt[ok, keep, drop = FALSE], Zt, rb, rep(0, n),
                     rep(1, n), TRUE, seq_along(keep), f2)),
    error = function(e) NULL)
  if (is.null(m) || !all(is.finite(m$theta))) return(NULL)
  m
}

#' The linear predictor a block's coefficients enter, by name, or `""`.
#'
#' @noRd
gp_block_lp <- function(bk, frame) {
  for (key in names(frame[["linpreds"]])) {
    Z <- frame[["linpreds"]][[key]][["Z"]]
    if (NCOL(Z) < max(bk[["c_idx"]]) || !NROW(Z)) next
    if (Matrix::nnzero(Z[, bk[["c_idx"]], drop = FALSE])) return(key)
  }
  ""
}

#' The response on the scale of its mean's link, clamped at the family's edge, or
#' `NULL` for a family with no such reading (ordinal, categorical, a mixture).
#'
#' @noRd
gp_working_y <- function(frame, lp) {
  key <- as.character(lp[["resp"]])
  y <- frame[["y"]][[key]]
  if (is.null(y) || !is.numeric(y) || !is.null(dim(y))) return(NULL)
  y <- as.numeric(y)
  fam <- tolower(frame[["spec"]]$responses[[key]]$family[["family"]] %||% "")
  link <- lp[["link"]]$linkfun
  if (fam %in% c("lognormal", "shifted_lognormal")) {
    return(gp_finite(log(y)))
  }
  m <- if (fam %in% c("poisson", "negbinomial", "negbinomial2", "geometric")) {
    y + 0.5
  } else if (fam == "binomial") {
    tr <- frame[["aterm_values"]][[key]][["trials"]]
    if (is.null(tr)) return(NULL)
    (y + 0.5) / (as.numeric(tr) + 1)
  } else if (fam == "bernoulli") {
    0.25 + 0.5 * y
  } else if (fam %in% c("gaussian", "student", "skew_normal", "exgaussian",
                        "asym_laplace", "huber", "gen_extreme_value", "gamma",
                        "weibull", "frechet", "exponential", "inverse.gaussian",
                        "beta")) {
    y
  } else {
    return(NULL)
  }
  gp_finite(suppressWarnings(link(m)))
}

#' A link's value with what fell off its domain at the nearest value it read. A
#' missing response stays missing.
#'
#' @noRd
gp_finite <- function(v) {
  ok <- is.finite(v)
  if (!any(ok)) return(NULL)
  v[!ok & !is.na(v) & v < 0] <- min(v[ok])
  v[!ok & !is.na(v) & v > 0] <- max(v[ok])
  v[is.nan(v)] <- min(v[ok])
  v
}

#' The outer start with each exact gp() at the marginal's maximum where the
#' marginal applies, at the working response's value and length scale where it
#' does not, and for what neither reaches at the best of a few length scales read
#' from the objective. `start` replaces what it names. A `theta` there turns all
#' of this off, and a `beta` or `betad` replaces only that component: the same
#' spelling with `start = list(beta = 5)` stopped 52 nats short when it threw the
#' warm start away.
#'
#' `to_beta` carries a natural `beta` to the coordinates the tape holds it in.
#' `body = FALSE` turns off the start of a block in a nonlinear body
#' (gp_body_start()), for a caller's own `start`. A restart from a saddle is
#' not the caller's, and keeps it.
#'
#' @noRd
gp_start_scan <- function(obj, template, frame, start, to_beta = identity,
                          body = TRUE) {
  bks <- Filter(function(bk) identical(bk[["covstruct"]], "gp"),
                frame[["re_blocks"]] %||% list())
  if (!length(bks) || !is.null(start[["theta"]])) return(obj)
  val <- function(tpl) {
    p <- outer_from_template(tpl, obj, frame)
    if (is.null(p)) return(NA_real_)
    v <- tryCatch(obj$fn(p), error = function(e) NA_real_)
    if (isTRUE(is.finite(v))) v else NA_real_
  }
  m <- gp_marginal_max(frame)
  if (!is.null(m)) {
    tpl <- template
    if (is.null(start[["beta"]])) tpl[["beta"]] <- to_beta(m$beta)
    tpl[["theta"]][] <- m$theta
    if (!is.null(m$betad) && is.null(start[["betad"]])) {
      tpl[["betad"]][] <- m$betad
    }
    if (!is.na(val(tpl))) {
      obj$par <- outer_from_template(tpl, obj, frame)
      # the proof check_convergence() reads: a stop at this value is at a maximum
      # the marginal found by another route
      obj$env$.frmtmb_gp_marginal <- m$value
      return(obj)
    }
  }
  w <- tryCatch(gp_start_working(frame, template, body = body),
                error = function(e) NULL)
  covered <- integer(0)
  scanned <- template
  if (!is.null(w)) {
    tpl <- template
    tpl[["theta"]][w$idx] <- w$theta
    if (!is.na(val(tpl))) {
      template <- tpl
      covered <- w$idx
    }
  }
  out <- gp_start_scan_units(obj, template, frame, val, covered)
  body <- intersect(covered, w[["body"]])
  if (!length(body)) return(out)
  # 0.79.1's start, the scan with the body's blocks not covered, is the second
  # descent's: on 7 of 192 rows it reached up to 3.30 nats more than the body
  # start (dev/reviews/2026-11-06-fitedge.md, B1). The inner solve starts where
  # the last evaluation left it, so reading it must not move the first descent
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  keep <- setdiff(covered, body)
  scanned[["theta"]][keep] <- template[["theta"]][keep]
  alt <- gp_start_scan_units(obj, scanned, frame, val, keep)$par
  if (length(alt) == length(out$par) && all(is.finite(alt))) {
    out$env$.frmtmb_gp_alt <- c(out$env$.frmtmb_gp_alt, list(alt))
  }
  out
}

#' A second descent from the working marginal's next local maximum
#' (gp_start_scan()), where it has one, and the better optimum kept. The working
#' response is not the model, and its best maximum can start the descent in a
#' basin of the model's that is not the best: on bernoulli, binomial and
#' low-count Poisson rows the fit stopped 0.2 to 1.07 nats short at conv 0, and
#' the next maximum reached the best of 15 starts (dev/arstep-findings.md).
#'
#' @noRd
gp_second_descent <- function(obj, opt, control, bounds, par_units = NULL,
                              verbose = 0L, tally = NULL) {
  alt <- obj$env$.frmtmb_gp_alt
  if (!length(alt)) return(opt)
  rm(".frmtmb_gp_alt", envir = obj$env)
  if (!all(is.finite(opt$par)) || !is.finite(opt$objective)) return(opt)
  saved <- obj_state_save(obj)
  best <- opt
  for (p in alt) {
    o <- tryCatch(optimize_obj(obj, control, bounds, par_units,
                               start_par = p, tally = tally),
                  error = function(e) NULL)
    if (!is.null(o) && all(is.finite(o$par)) && is.finite(o$objective) &&
          best$objective - o$objective > 1e-10 * max(1, abs(best$objective))) {
      best <- o
    }
  }
  obj_settle(obj, opt, best, saved)
  gain <- opt$objective - best$objective
  if (verbose) {
    vb_say("exact gp(): ", length(alt), " more start(s), gained ",
           format(gain, digits = 3), " in log-likelihood")
  }
  best$gp_second <- c(runs = length(alt), gain = gain)
  # every optimizer run is in the count, the lost one too: se_check_at_fit()
  # weighs the fit-time Hessian against it
  if (!is.null(tally)) {
    best$evals <- tally$evals
    if (!is.null(best$nonfinite_trials)) best$nonfinite_trials <- tally$n
  }
  best
}

#' The best of a few length scales, read from the objective with each field's
#' sd on its predictor's unit and, where the response has a sigma, sigma on a
#' ladder below that. In y's units, because the objective at unit sd against a y
#' of size 1000 reads 4e6 on every length scale and tells them apart by nothing.
#'
#' A field that is rough against its positions has a second maximum at a long
#' length scale with the field near zero, which a descent from a long start ends
#' in. On a predictor whose link has no unit, the unit is 1 and there is no ladder.
#'
#' @noRd
gp_start_scan_units <- function(obj, template, frame, val, covered = integer(0)) {
  bks <- Filter(function(bk) {
    identical(bk[["covstruct"]], "gp") && !any(bk[["theta_idx"]] %in% covered)
  }, frame[["re_blocks"]] %||% list())
  units <- tryCatch(predictor_units(frame), error = function(e) numeric(0))
  best <- template
  fb <- val(best)
  if (is.na(fb)) fb <- Inf
  for (bk in bks) {
    pos <- bk[["theta_idx"]][-1L]
    d <- sqrt(max(Reduce(`+`, bk[["aux_D2"]])))
    if (!isTRUE(d > 0)) next
    key <- gp_block_lp(bk, frame)
    u <- if (key %in% names(units)) units[[key]] else 1
    sg <- gp_scan_sigma_slot(frame, bk)
    ladder <- if (is.null(sg)) NA_real_ else u * c(1, 0.1, 1e-2, 1e-3)
    sg_link <- attr(sg, "linkfun")
    for (r in d / 4^(0:3)) {
      for (nv in ladder) {
        tpl <- best
        tpl[["theta"]][pos] <- log(r)
        tpl[["theta"]][bk[["theta_idx"]][1L]] <- log(u)
        if (!is.na(nv)) tpl[["betad"]][sg] <- sg_link(nv)
        v <- val(tpl)
        if (!is.na(v) && v < fb) {
          fb <- v
          best <- tpl
        }
      }
    }
  }
  p <- outer_from_template(best, obj, frame)
  if (!is.null(p)) obj$par <- p
  obj
}

#' Where a block's response keeps a modeled sigma's intercept in `betad`, with
#' its link as an attribute, or `NULL`. On a block of sigma itself, or on a family
#' without one, there is nothing to ladder.
#'
#' @noRd
gp_scan_sigma_slot <- function(frame, bk) {
  lp <- frame[["linpreds"]][[gp_block_lp(bk, frame)]]
  if (is.null(lp) || !identical(lp[["dpar"]], "mu")) return(NULL)
  for (l in frame[["linpreds"]]) {
    if (!identical(l[["dpar"]], "sigma") || !identical(l[["par"]], "betad") ||
          !identical(l[["resp"]], lp[["resp"]])) {
      next
    }
    icpt <- match("(Intercept)", colnames(l[["X"]]))
    if (is.na(icpt) || !is.null(l[["constant"]])) return(NULL)
    if (l[["idx"]][icpt] %in% frame[["betad_fixed_idx"]]) return(NULL)
    return(structure(l[["idx"]][icpt], linkfun = l[["link"]]$linkfun))
  }
  NULL
}
