# Lane `latent`: the two slicing rules and the profile arithmetic, shared
# by dev/latent-2p4-startfix.R and dev/latent-2p4-shrink.R so that the
# weight sweep measures the SAME construction the comparison did.

# today's: the mean item code, rank-cut into K equal-count slices
slice_score <- function(y, ncat, K) {
  n <- nrow(y)
  span <- pmax(ncat - 1L, 1L)
  sc <- rowMeans(sweep(y - 1, 2L, span, "/"), na.rm = TRUE)
  fin <- is.finite(sc)
  sc[!fin] <- if (any(fin)) stats::median(sc[fin]) else 0
  pmin(1L + ((rank(sc, ties.method = "first") - 1L) * K) %/% n, K)
}

# the candidate: hard k-means on the one-hot response pattern, seeded by
# farthest point, relabelled by the same mean score
slice_kmeans <- function(y, ncat, K, iter = 25L) {
  n <- nrow(y)
  J <- ncol(y)
  X <- matrix(0, n, sum(ncat))
  off <- 0L
  for (j in seq_len(J)) {
    v <- y[, j]
    md <- which.max(tabulate(v[!is.na(v)], nbins = ncat[j]))
    v[is.na(v)] <- md
    X[cbind(seq_len(n), off + v)] <- 1
    off <- off + ncat[j]
  }
  mu0 <- colMeans(X)
  d2 <- rowSums(sweep(X, 2L, mu0, "-")^2)
  idx <- integer(K)
  idx[1L] <- which.max(d2)
  if (K > 1L) {
    dmin <- rowSums(sweep(X, 2L, X[idx[1L], ], "-")^2)
    for (k in 2:K) {
      idx[k] <- which.max(dmin)
      dmin <- pmin(dmin, rowSums(sweep(X, 2L, X[idx[k], ], "-")^2))
    }
  }
  C <- X[idx, , drop = FALSE]
  cl <- rep(1L, n)
  for (it in seq_len(iter)) {
    dd <- vapply(seq_len(K),
                 function(k) rowSums(sweep(X, 2L, C[k, ], "-")^2),
                 numeric(n))
    new <- max.col(-dd, ties.method = "first")
    if (identical(new, cl) && it > 1L) break
    cl <- new
    for (k in seq_len(K)) {
      if (any(cl == k)) C[k, ] <- colMeans(X[cl == k, , drop = FALSE])
    }
  }
  # a cluster that emptied keeps its seed row, so every class starts
  # from a real response pattern rather than from a degenerate profile
  for (k in seq_len(K)) if (!any(cl == k)) cl[idx[k]] <- k
  # relabel by the same mean score the shipped rule uses, so the
  # documented low-to-high class ordering survives where it means
  # anything and is at least deterministic where it does not
  span <- pmax(ncat - 1L, 1L)
  sc <- rowMeans(sweep(y - 1, 2L, span, "/"), na.rm = TRUE)
  sc[!is.finite(sc)] <- 0
  ord <- order(vapply(seq_len(K), function(k) mean(sc[cl == k]), 0))
  match(cl, ord)
}

# the profile logits a slicing produces, the shipped arithmetic exactly,
# with an optional shrink toward the POOLED profile. `w = 0` is the
# shipped rule; a positive `w` keeps the classes separated but stops the
# start from being as confident as a hard partition makes it.
profiles_from <- function(y, ncat, K, slice, w = 0) {
  J <- ncol(y)
  out <- vector("list", J)
  for (j in seq_len(J)) {
    aj <- y[, j]
    aj <- aj[!is.na(aj)]
    pool <- (tabulate(aj, nbins = ncat[j]) + 1) /
      sum(tabulate(aj, nbins = ncat[j]) + 1)
    v <- numeric(0)
    for (k in seq_len(K)) {
      cj <- y[slice == k, j]
      cj <- cj[!is.na(cj)]
      cnt <- tabulate(cj, nbins = ncat[j]) + 1
      p <- cnt / sum(cnt)
      p <- (1 - w) * p + w * pool
      v <- c(v, log(p[-1L]) - log(p[1L]))
    }
    out[[j]] <- v
  }
  stats::setNames(out, paste0("pi", seq_len(J)))
}

