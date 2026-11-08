# The joint covariance of a fit with a latent residual, held sparse.
#
# WHY. A latent residual puts an inner parameter on every row, so the dense
# joint covariance `V = Q^-1` of get_joint_cov() is `(n + p)^2` numbers: 32 MB
# kept on the fit at 2000 rows, 0.5 GB at 8000, 8.2 GB at 32000. Its product
# with the dense design of the fitted rows cost 4.7 s at 2000 rows and 332 s
# at 8000 (dev/latperf-findings.md). The joint precision `Q` is sparse, banded
# along an AR chain, and no reader needs `V` whole: a new row reads a block of
# `V` at a few positions (`beta`, `theta` and the levels of another block), and
# the fitted rows read `diag(A V A')`, with `A` sparse.
#
# THE ALGEBRA. With `r` the inner positions and `f` the outer ones,
#
#   V_rr = W + J Vf J',  V_rf = J Vf,  V_ff = Vf,
#   W = Q_rr^-1,  J = -W Q_rf,
#
# which is the block inverse joint_cov_repair() writes. `Vf` is the covariance
# the standard-error check kept for propagation where it kept one, and
# otherwise the inverse of the Schur complement `Q_ff + Q_fr J`, which is the
# `f` block of `solve(Q)`. A row `a = (a_r, a_f)` then has
#
#   a' V a = a_r' W a_r + g' Vf g,   g = J' a_r + a_f.
#
# `a_r' W a_r` reads `W` only at pairs of positions that one row loads
# together. Where the row's likelihood couples a pair in `Q_rr`, the pair is in
# the pattern of the Cholesky factor of `Q_rr`, where the selected inverse is
# exact: TMB's, the one sdreport() uses for `diag.cov.random`. A family that
# separates its dpars does not couple every such pair: a hurdle's mean
# `(1 - hu) mu` loads the latent of `mu` and a level of `hu`, which no row's
# likelihood joins. Those pairs are read from columns of `W` solved for once
# per distinct position on the shorter side, a level of `hu` here, so the cost
# stays linear in the rows (jcs_W_pairs()).
#
# Everything kept is O(n): the factor, `J` (n x |f|) and `Vf`. The selected
# inverse is rebuilt by each reader that needs it, 0.14 s at 32000 rows, so a
# saved fit does not carry it.

#' The most inner positions a dense block of `W` is formed for by direct
#' solves, one column each. A new row loads the levels of another block, a
#' few dozen; a fitted row loads every latent, and takes the selected inverse.
#'
#' @noRd
jcs_block_max <- 200L

#' Whether a joint covariance (get_joint_cov()) is held sparse.
#'
#' @noRd
jc_is_sparse <- function(jc) inherits(jc[["V"]], "frm_joint_sparse")

#' The joint covariance of a fit with a latent residual, held sparse, in the
#' shape get_joint_cov() returns, or `NULL` where this route does not apply:
#' a covariance held as a factor (se_reparam_cov(), sdr_unwhiten()), an inner
#' Hessian that does not factor, or a Schur complement that does not invert.
#' get_joint_cov() then takes the dense route, which says why.
#'
#' @noRd
jcs_joint_cov <- function(fit) {
  Q <- joint_precision(fit)
  if (is.null(Q)) return(NULL)
  sdr <- sdr_of(fit)
  if (!is.null(sdr$cov_fixed_root) || !is.null(sdr$joint_root)) return(NULL)
  p <- nrow(Q)
  r <- fit$obj$env$random
  f <- setdiff(seq_len(p), r)
  if (!length(r) || !length(f)) return(NULL)
  rn <- rownames(Q)
  Q <- methods::as(methods::as(Q, "CsparseMatrix"), "generalMatrix")
  Qrr <- Matrix::forceSymmetric(Q[r, r, drop = FALSE])
  # CHOLMOD narrates a failed factorization from C; the dense route below
  # gives the message that names the remedy
  L <- tryCatch(suppressWarnings(Matrix::Cholesky(Qrr, perm = TRUE,
                                                  super = TRUE, LDL = FALSE)),
                error = function(e) NULL)
  if (is.null(L)) return(NULL)
  Qrf <- Q[r, f, drop = FALSE]
  J <- -as.matrix(Matrix::solve(L, Qrf))
  if (!all(is.finite(J))) return(NULL)
  Vp <- sdr$cov_fixed_prop
  # joint_cov_repair()'s conditions, so both routes read the same `Vf`
  repair <- !is.null(Vp) && length(f) == nrow(Vp) &&
    length(f) == length(outer_par_names(fit))
  if (repair) {
    Vf <- unname(as.matrix(Vp))
    meta <- joint_cov_repair_meta(fit, sdr, p, r, f)
  } else {
    S <- as.matrix(Q[f, f, drop = FALSE]) +
      as.matrix(Matrix::crossprod(Qrf, J))
    # base solve() refuses below its reciprocal-condition limit, and the fit
    # then takes the dense route, which inverts the whole `Q` as before: an
    # `unstr()` fit with 95 lost standard errors reads rcond 1.6e-101 here
    Vf <- tryCatch(solve(S), error = function(e) NULL)
    if (is.null(Vf) || !all(is.finite(Vf))) return(NULL)
    meta <- list(null = NULL, flat = NULL, units = NULL, lost_pos = integer(0))
  }
  V <- structure(list(L = L, r = r, f = f, J = unname(J), Vf = unname(Vf),
                      p = p, names = rn),
                 class = "frm_joint_sparse")
  lost <- meta$lost_pos
  list(V = V, names = rn, null = meta$null, flat = meta$flat,
       units = meta$units, lost_pos = lost[!is.na(lost)], root = NULL,
       Vrest = NULL)
}

# A block of a sparse joint covariance as a dense matrix: `V[i, j]` at
# positions `i` and `j`, with the dimnames the dense route gives. Every reader
# of a block of `jc$V` reaches this, so a block of a few positions costs a few
# solves. A block of more than `frmtmb.latres_dense_rows`^2 numbers is refused,
# because it is the dense matrix this route avoids.
#' @export
`[.frm_joint_sparse` <- function(x, i, j, ..., drop = TRUE) {
  frm_check_dots(...)
  p <- x[["p"]]
  i <- if (missing(i)) seq_len(p) else seq_len(p)[i]
  j <- if (missing(j)) seq_len(p) else seq_len(p)[j]
  jcs_check_size(length(i), length(j))
  r <- x[["r"]]
  f <- x[["f"]]
  ri <- match(i, r)
  rj <- match(j, r)
  fi <- match(i, f)
  fj <- match(j, f)
  a <- which(!is.na(ri))
  b <- which(!is.na(rj))
  af <- which(!is.na(fi))
  bf <- which(!is.na(fj))
  nm <- x[["names"]]
  out <- matrix(0, length(i), length(j),
                dimnames = if (!is.null(nm)) list(nm[i], nm[j]))
  JV <- x[["J"]][ri[a], , drop = FALSE] %*% x[["Vf"]]
  if (length(a) && length(b)) {
    out[a, b] <- jcs_W(x, ri[a], rj[b]) +
      tcrossprod(JV, x[["J"]][rj[b], , drop = FALSE])
  }
  if (length(a) && length(bf)) out[a, bf] <- JV[, fj[bf], drop = FALSE]
  if (length(af) && length(b)) {
    out[af, b] <- t(x[["J"]][rj[b], , drop = FALSE] %*%
                      x[["Vf"]][, fi[af], drop = FALSE])
  }
  if (length(af) && length(bf)) {
    out[af, bf] <- x[["Vf"]][fi[af], fj[bf], drop = FALSE]
  }
  if (drop) drop(out) else out
}

# The whole matrix, for frm_joint_cov(), under the same limit.
#' @export
as.matrix.frm_joint_sparse <- function(x, ...) {
  frm_check_dots(...)
  x[seq_len(x[["p"]]), seq_len(x[["p"]]), drop = FALSE]
}

#' The refusal of a dense block of `ni x nj` numbers above the limit.
#'
#' @noRd
jcs_check_size <- function(ni, nj) {
  lim <- getOption("frmtmb.latres_dense_rows", 20000L)
  if (as.numeric(ni) * nj <= as.numeric(lim)^2) return(invisible(NULL))
  gb <- 8 * as.numeric(ni) * nj / 1e9
  # the solved columns, J Vf J', the block and its copy: the peak was 4.0 to
  # 4.1 times the block at 1204 and 3004 positions
  frm_stop("This reads a block of the joint covariance of a latent residual ",
           "as a dense matrix of ", ni, " x ", nj, ", about ",
           format(gb, digits = 2), " GB, and about ", format(4 * gb, digits = 2),
           " GB while it is formed. ",
           "options(frmtmb.latres_dense_rows = ) moves the limit of ", lim,
           " rows; frm_linpred(), predict(), and fitted() and residuals() ",
           "of a family with one mean per row do not read it", call. = FALSE,
           class = "frmtmb_latres_dense")
}

#' `W[ri, rj]`, `W = Q_rr^-1`, dense, by one solve per distinct position on
#' the shorter side.
#'
#' @noRd
jcs_W <- function(S, ri, rj) {
  if (length(unique(rj)) > length(unique(ri))) return(t(jcs_W(S, rj, ri)))
  u <- unique(rj)
  E <- Matrix::sparseMatrix(i = u, j = seq_along(u), x = 1,
                            dims = c(length(S[["r"]]), length(u)))
  X <- as.matrix(Matrix::solve(S[["L"]], E))
  X[ri, match(rj, u), drop = FALSE]
}

#' The selected inverse of `Q_rr` on the pattern of its supernodal factor, or
#' `NULL` where TMB does not provide it.
#'
#' @noRd
jcs_selinv <- function(S) {
  # TMB's own routine (C, the one sdreport() reads diag.cov.random with);
  # unexported, so reached by name and given up on if it is gone
  ss <- get0("solveSubset", envir = asNamespace("TMB"), inherits = FALSE)
  if (!is.function(ss)) return(NULL)
  tryCatch(ss(L = S[["L"]]), error = function(e) NULL)
}

#' `A` as a column-compressed sparse matrix with no column that no row loads
#' dropped, and the positions that go with its columns.
#'
#' @noRd
jcs_trim <- function(A, pos) {
  A <- methods::as(methods::as(methods::as(A, "CsparseMatrix"),
                               "generalMatrix"), "dMatrix")
  keep <- diff(A@p) > 0L
  list(A = A[, keep, drop = FALSE], pos = pos[keep])
}

#' Which columns of `A`, dense or sparse, any row loads: a stored entry, or
#' a value that is not exactly zero.
#'
#' @noRd
jcs_loaded_cols <- function(A) {
  if (inherits(A, "Matrix")) {
    A <- methods::as(methods::as(A, "CsparseMatrix"), "generalMatrix")
    return(diff(A@p) > 0L)
  }
  colSums(A != 0 | is.na(A)) > 0
}

#' `g = J' a_r + a_f` for every row of `A` over positions `pos`: the rows'
#' load on the outer parameters, through the inner ones.
#'
#' @noRd
jcs_outer_load <- function(S, A, pos) {
  ri <- match(pos, S[["r"]])
  fi <- match(pos, S[["f"]])
  kr <- which(!is.na(ri))
  G <- matrix(0, nrow(A), length(S[["f"]]))
  if (length(kr)) {
    G <- G + as.matrix(A[, kr, drop = FALSE] %*%
                         S[["J"]][ri[kr], , drop = FALSE])
  }
  for (k in which(!is.na(fi))) {
    G[, fi[k]] <- G[, fi[k]] + as.numeric(A[, k])
  }
  G
}

#' `diag(A V[pos, pos] A')` of a sparse joint covariance `S`.
#'
#' @noRd
jcs_rowvar <- function(S, A, pos) {
  tr <- jcs_trim(A, pos)
  A <- tr$A
  pos <- tr$pos
  G <- jcs_outer_load(S, A, pos)
  v <- rowSums((G %*% S[["Vf"]]) * G)
  ri <- match(pos, S[["r"]])
  kr <- which(!is.na(ri))
  if (length(kr)) {
    v <- v + jcs_rowvar_inner(S, A[, kr, drop = FALSE], ri[kr])
  }
  v
}

#' `diag(Ar W[ri, ri] Ar')`: by a dense block of `W` where the rows load few
#' inner positions, by direct solves where there are few rows, and otherwise
#' through the selected inverse, with each pair outside its pattern read from
#' a column of `W` solved once (jcs_W_pairs()).
#'
#' @noRd
jcs_rowvar_inner <- function(S, Ar, ri) {
  if (length(unique(ri)) <= jcs_block_max) {
    Wb <- jcs_W(S, ri, ri)
    return(as.numeric(Matrix::rowSums((Ar %*% Wb) * Ar)))
  }
  if (nrow(Ar) <= jcs_block_max) return(jcs_rowvar_solve(S, Ar, ri))
  W <- jcs_selinv_keys(S)
  if (is.null(W)) return(jcs_rowvar_solve(S, Ar, ri))
  # one term per pair of a row's own entries, looked up in the pattern:
  # `(Ar %*% W) * Ar` would be dense, since a level another block shares
  # couples its column of W to every latent
  t <- methods::as(Ar, "TsparseMatrix")
  o <- order(t@i)
  row <- t@i[o] + 1L
  k <- ri[t@j[o] + 1L]
  x <- t@x[o]
  m <- tabulate(row, nbins = nrow(Ar))
  mi <- m[row]
  e1 <- rep(seq_along(row), times = mi)
  e2 <- sequence(mi, from = cumsum(c(1L, m))[row])
  a <- pmin(k[e1], k[e2]) - 1
  b <- pmax(k[e1], k[e2]) - 1
  hit <- match(a * W$N + b, W$key)
  w <- W$x[hit]
  # a pair outside the pattern is a structural zero where W is not zero
  out <- which(is.na(hit))
  if (length(out)) w[out] <- jcs_W_pairs(S, a[out] + 1, b[out] + 1)
  term <- x[e1] * x[e2] * w
  v <- numeric(nrow(Ar))
  s <- rowsum(term, row[e1], reorder = TRUE)
  v[as.integer(rownames(s))] <- s[, 1L]
  v
}

#' `W[i[k], j[k]]` for each pair `k`, from columns of `W` solved once each.
#' Each pair takes its column at the end that more pairs share: a hurdle's
#' latent-and-level pairs take the level's, so 50 levels cost 50 solves
#' whatever the rows. The columns are solved in chunks of about 32 MB. The
#' attribute `columns` is how many were solved.
#'
#' @noRd
jcs_W_pairs <- function(S, i, j) {
  nr <- length(S[["r"]])
  freq <- tabulate(c(i, j), nbins = nr)
  by_i <- freq[i] >= freq[j]
  col <- ifelse(by_i, i, j)
  oth <- ifelse(by_i, j, i)
  u <- unique(col)
  out <- numeric(length(i))
  size <- max(1L, floor(2^22 / nr))
  for (blk in split(u, ceiling(seq_along(u) / size))) {
    E <- Matrix::sparseMatrix(i = blk, j = seq_along(blk), x = 1,
                              dims = c(nr, length(blk)))
    X <- as.matrix(Matrix::solve(S[["L"]], E))
    here <- which(col %in% blk)
    out[here] <- X[cbind(oth[here], match(col[here], blk))]
  }
  attr(out, "columns") <- length(u)
  out
}

#' The selected inverse as a lookup: the key `min * N + max` of each stored
#' entry (0-based) and its value. Rebuilt per reader, not memoized: it is the
#' largest part of the sparse covariance, and 0.14 s at 32000 rows.
#'
#' @noRd
jcs_selinv_keys <- function(S) {
  W <- jcs_selinv(S)
  if (is.null(W)) return(NULL)
  t <- methods::as(W, "TsparseMatrix")
  N <- as.numeric(nrow(W))
  list(N = N, key = pmin(t@i, t@j) * N + pmax(t@i, t@j), x = t@x)
}


#' `diag(Ar W[ri, ri] Ar')` by solving `Q_rr x = a_r` for each row, in
#' chunks of about 32 MB.
#'
#' @noRd
jcs_rowvar_solve <- function(S, Ar, ri) {
  nr <- length(S[["r"]])
  E <- Matrix::sparseMatrix(i = ri, j = seq_along(ri), x = 1,
                            dims = c(nr, length(ri)))
  n <- nrow(Ar)
  out <- numeric(n)
  size <- max(1L, floor(2^22 / nr))
  for (blk in split(seq_len(n), ceiling(seq_len(n) / size))) {
    B <- E %*% Matrix::t(Ar[blk, , drop = FALSE])
    X <- as.matrix(Matrix::solve(S[["L"]], B))
    out[blk] <- as.numeric(Matrix::colSums(B * X))
  }
  out
}

#' `G V[pos, pos] G'` of a sparse joint covariance, for a few rows of `G`.
#'
#' @noRd
jcs_quad <- function(S, G, pos) {
  tr <- jcs_trim(G, pos)
  A <- tr$A
  pos <- tr$pos
  Gf <- jcs_outer_load(S, A, pos)
  out <- Gf %*% S[["Vf"]] %*% t(Gf)
  ri <- match(pos, S[["r"]])
  kr <- which(!is.na(ri))
  if (length(kr)) {
    Ar <- A[, kr, drop = FALSE]
    rk <- ri[kr]
    if (length(unique(rk)) <= jcs_block_max) {
      out <- out + as.matrix(Ar %*% jcs_W(S, rk, rk) %*% Matrix::t(Ar))
    } else {
      E <- Matrix::sparseMatrix(i = rk, j = seq_along(rk), x = 1,
                                dims = c(length(S[["r"]]), length(rk)))
      B <- E %*% Matrix::t(Ar)
      out <- out + as.matrix(Matrix::crossprod(B, Matrix::solve(S[["L"]], B)))
    }
  }
  out
}

#' predict_mean_se()'s gradient over the union of coefficient positions
#' `pos`, sparse: each dpar's design scaled by its row's `d mean / d eta`,
#' added where two dpars share a position. Dense, over every latent, it is
#' `n^2` numbers.
#'
#' @noRd
jcs_stack_grad <- function(das, grad, dnames, pos, n) {
  G <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                            dims = c(n, length(pos)))
  for (dnm in dnames) {
    cp <- das[[dnm]]$coef_pos
    if (!length(cp)) next
    M <- Matrix::sparseMatrix(i = seq_along(cp), j = match(cp, pos), x = 1,
                              dims = c(length(cp), length(pos)))
    A <- methods::as(das[[dnm]]$A, "CsparseMatrix")
    G <- G + Matrix::Diagonal(x = grad[[dnm]]) %*% A %*% M
  }
  G
}

#' `diag(V)[idx]` of a joint covariance, dense or sparse.
#'
#' @noRd
jc_diag <- function(jc, idx) {
  if (!jc_is_sparse(jc)) return(diag(jc$V)[idx])
  jcs_rowvar(jc$V, Matrix::Diagonal(length(idx)), idx)
}
