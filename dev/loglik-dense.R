# Lane loglik: the dense reference for frmtmb's leave-one-out log
# density. The joint covariance of the response is built from the
# fitted model's own definition (autocor_matrix() over the time levels
# for the ar() family, sigma^2 V for fcor(), sigma^2 (A'A)^-1 for sar()),
# and each row's conditional is taken from it by the textbook partition
# formula, with no precision matrix. Sourced by dev/loglik-*.R; reads
# `d` (with columns g and t for the ar() family), `V` and `W` from the
# caller.

cond_one <- function(y, m, S, nu = NULL) {
  n <- length(y)
  vapply(seq_len(n), function(i) {
    if (n == 1L) {
      return(if (is.null(nu)) dnorm(y, m, sqrt(S[1, 1]), log = TRUE) else
        dt((y - m) / sqrt(S[1, 1]), nu, log = TRUE) - 0.5 * log(S[1, 1]))
    }
    o <- seq_len(n)[-i]
    A <- solve(S[o, o, drop = FALSE])
    r <- y[o] - m[o]
    cm <- m[i] + as.numeric(S[i, o, drop = FALSE] %*% A %*% r)
    cv <- S[i, i] - as.numeric(S[i, o, drop = FALSE] %*% A %*% S[o, i])
    if (is.null(nu)) return(dnorm(y[i], cm, sqrt(cv), log = TRUE))
    d2 <- as.numeric(t(r) %*% A %*% r)
    df <- nu + n - 1
    sc <- sqrt((nu + d2) / df * cv)
    dt((y[i] - cm) / sc, df, log = TRUE) - log(sc)
  }, 0)
}

# The model's joint mean and covariance and its groups, from the fit.
dense_model <- function(fit, dd = d) {
  ac <- fit$frame$autocor[[1]]
  dp <- eval_dpars(fit)[[1]]
  y <- fit$frame$y[[1]]
  nn <- length(y)
  mu <- rep_len(as.numeric(dp$mu), nn)
  s <- rep_len(as.numeric(dp$sigma), nn)
  nu <- if (!is.null(dp$nu)) as.numeric(dp$nu)[1]
  if (identical(ac$struct, "fcor")) {
    return(list(y = y, m = mu, S = outer(s, s) * V, grp = rep(1, nn),
                nu = nu))
  }
  if (identical(ac$struct, "sar")) {
    rho <- frmtmb:::autocor_natural(fit$estimates$thetaac[ac$theta_idx], ac)
    A <- diag(nn) - rho * W
    m <- if (ac$aux$type == "lag") as.vector(solve(A, mu)) else mu
    return(list(y = y, m = m, S = s[1]^2 * solve(crossprod(A)),
                grp = rep(1, nn), nu = nu))
  }
  R <- autocor_matrix(fit)
  ti <- match(as.character(dd$t), rownames(R))
  S <- matrix(0, nn, nn)
  for (g in unique(dd$g)) {
    r <- which(dd$g == g)
    S[r, r] <- outer(s[r], s[r]) * R[ti[r], ti[r]]
  }
  list(y = y, m = mu, S = S, grp = dd$g, nu = nu)
}

dense_loo <- function(fit, dd = d) {
  dm <- dense_model(fit, dd)
  out <- numeric(length(dm$y))
  for (g in unique(dm$grp)) {
    r <- which(dm$grp == g)
    out[r] <- cond_one(dm$y[r], dm$m[r], dm$S[r, r, drop = FALSE], dm$nu)
  }
  out
}

# The joint log density of the same model, for the check that the row
# sum is not it.
dense_joint <- function(fit, dd = d) {
  dm <- dense_model(fit, dd)
  tot <- 0
  for (g in unique(dm$grp)) {
    r <- which(dm$grp == g)
    S <- dm$S[r, r, drop = FALSE]
    e <- dm$y[r] - dm$m[r]
    k <- length(r)
    q <- sum(e * solve(S, e))
    ld <- as.numeric(determinant(S)$modulus)
    tot <- tot + if (is.null(dm$nu)) {
      -0.5 * k * log(2 * pi) - 0.5 * ld - 0.5 * q
    } else {
      nu <- dm$nu
      lgamma((nu + k) / 2) - lgamma(nu / 2) - 0.5 * k * log(nu * pi) -
        0.5 * ld - 0.5 * (nu + k) * log1p(q / nu)
    }
  }
  tot
}
