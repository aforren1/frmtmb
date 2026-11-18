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

#' Whether a REML fit on a whitened design (beta_whitener()) reads its
#' joint covariance by sparse solves in the whitened coordinates: above
#' joint_root_limit() joint parameters, or where sdr_unwhiten() kept some
#' rows of the factor only (joint_root_rows()).
#'
#' @noRd
jcs_white_wanted <- function(fit) {
  if (is.null(fit$beta_whiten)) return(FALSE)
  sdr <- tryCatch(sdr_of(fit), error = keep_missing(NULL))
  Qg <- sdr$joint_precision_white
  !is.null(Qg) && (nrow(Qg) > joint_root_limit() ||
                     length(attr(sdr$joint_root, "rows")) > 0L)
}

#' The joint covariance of a fit with a latent residual, held sparse, in the
#' shape get_joint_cov() returns, or `NULL` where this route does not apply:
#' a covariance held as a factor (se_reparam_cov(), sdr_unwhiten()), an inner
#' Hessian that does not factor, or a Schur complement that does not invert.
#' get_joint_cov() then takes the dense route, which says why.
#'
#' With `white`, the covariance of a REML fit on a whitened design
#' (jcs_white_wanted()), held in the whitened coordinates of the precision
#' sdr_unwhiten() was handed, `natural = M gamma`, with `M` kept
#' (joint_whitener()): a reader carries its loads there first,
#' `a' V a = (M' a)' V_gamma (M' a)` (jcs_white_load()). In the natural
#' coordinates `beta`'s block has the design's condition number squared
#' (joint_cov_white_inner()). There `Vf` keeps its factor where it has one
#' (`Rf`, se_reparam_cov()). `whole` takes the inverse of the whole
#' precision even where the standard-error check repaired `Vf`, as
#' joint_cov_whole() does.
#'
#' @noRd
jcs_joint_cov <- function(fit, white = jcs_white_wanted(fit), whole = FALSE) {
  sdr <- sdr_of(fit)
  Q <- if (white) sdr$joint_precision_white else joint_precision(fit)
  if (is.null(Q)) return(NULL)
  if (!white && (!is.null(sdr$cov_fixed_root) || !is.null(sdr$joint_root))) {
    return(NULL)
  }
  p <- nrow(Q)
  r <- fit$obj$env$random
  f <- setdiff(seq_len(p), r)
  if (!length(r) || !length(f)) return(NULL)
  rn <- rownames(Q)
  M <- NULL
  if (white) {
    jb <- which(rn == "beta")
    if (length(jb) != length(fit$frame[["par_template"]][["beta"]]) ||
          anyNA(match(jb, r))) {
      return(NULL)
    }
    M <- joint_whitener(fit$beta_whiten, jb, p)
  }
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
  repair <- !whole && !is.null(Vp) && length(f) == nrow(Vp) &&
    length(f) == length(outer_par_names(fit))
  Rf <- NULL
  if (repair) {
    Vf <- unname(as.matrix(Vp))
    meta <- joint_cov_repair_meta(fit, sdr, p, r, f)
    Lf <- sdr$cov_fixed_root
    if (white && !is.null(Lf) && nrow(Lf) == length(f)) Rf <- unname(Lf)
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
                      p = p, names = rn, M = M, Rf = Rf),
                 class = "frm_joint_sparse")
  lost <- meta$lost_pos
  why <- meta$lost_why
  list(V = V, names = rn, null = meta$null, flat = meta$flat,
       units = meta$units, lost_pos = lost[!is.na(lost)], root = NULL,
       Vrest = NULL,
       lost_why = if (length(why) == length(lost)) why[!is.na(lost)],
       limit_range = meta$limit_range)
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
  jcs_check_size(length(i), length(j), x)
  M <- x[["M"]]
  out <- if (is.null(M)) jcs_block(x, i, j) else {
    # held in the whitened coordinates: M[i, ] V_gamma M[j, ]'
    Mi <- M[i, , drop = FALSE]
    Mj <- M[j, , drop = FALSE]
    ci <- which(diff(Mi@p) > 0L)
    cj <- which(diff(Mj@p) > 0L)
    as.matrix(Mi[, ci, drop = FALSE] %*% jcs_block(x, ci, cj) %*%
                Matrix::t(Mj[, cj, drop = FALSE]))
  }
  nm <- x[["names"]]
  dimnames(out) <- if (!is.null(nm)) list(nm[i], nm[j])
  if (drop) drop(out) else out
}

#' `V[i, j]` of a sparse joint covariance in the coordinates it is held in.
#'
#' @noRd
jcs_block <- function(x, i, j) {
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
  out <- matrix(0, length(i), length(j))
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
  out
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
jcs_check_size <- function(ni, nj, x = NULL) {
  lim <- getOption("frmtmb.latres_dense_rows", 20000L)
  if (as.numeric(ni) * nj <= as.numeric(lim)^2) return(invisible(NULL))
  gb <- 8 * as.numeric(ni) * nj / 1e9
  # the solved columns, J Vf J', the block and its copy: the peak was 4.0 to
  # 4.1 times the block at 1204 and 3004 positions
  frm_stop("This reads a block of the joint covariance of ",
           if (is.null(x[["M"]])) "a latent residual" else "this REML fit",
           " as a dense matrix of ", ni, " x ", nj, ", about ",
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

#' Loads `A` over positions `pos` carried to the coordinates `S` is held in:
#' `A M[pos, ]` over every position, where `S` is held in whitened ones
#' (jcs_joint_cov()), else unchanged. Formed in the natural coordinates,
#' `a' V a` of a fitted value of a degree-7 raw polynomial cancels the
#' digits the whitener kept (joint_cov_white_inner()).
#'
#' @noRd
jcs_white_load <- function(S, A, pos) {
  M <- S[["M"]]
  if (is.null(M)) return(list(A = A, pos = pos))
  A <- methods::as(methods::as(methods::as(A, "CsparseMatrix"),
                               "generalMatrix"), "dMatrix")
  list(A = A %*% M[pos, , drop = FALSE], pos = seq_len(ncol(M)))
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
  wl <- jcs_white_load(S, A, pos)
  tr <- jcs_trim(wl$A, wl$pos)
  A <- tr$A
  pos <- tr$pos
  G <- jcs_outer_load(S, A, pos)
  v <- if (is.null(S[["Rf"]])) rowSums((G %*% S[["Vf"]]) * G) else {
    rowSums((G %*% S[["Rf"]])^2)
  }
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
  wl <- jcs_white_load(S, G, pos)
  tr <- jcs_trim(wl$A, wl$pos)
  A <- tr$A
  pos <- tr$pos
  Gf <- jcs_outer_load(S, A, pos)
  out <- if (is.null(S[["Rf"]])) Gf %*% S[["Vf"]] %*% t(Gf) else {
    tcrossprod(Gf %*% S[["Rf"]])
  }
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

#' A factor `L` of a covariance `V`, `V = L L'`, with a column per
#' direction it carries: a position with zero variance (a lost parameter
#' of a repaired covariance) gets a zero row, and a singular remainder
#' keeps the eigen-directions above 1e-12 of the largest, in the scale of
#' `V`'s own standard deviations.
#'
#' @noRd
jcs_psd_factor <- function(V) {
  d <- diag(V)
  on <- which(d > 0)
  L <- matrix(0, nrow(V), length(on))
  if (!length(on)) return(L)
  s <- sqrt(d[on])
  C <- V[on, on, drop = FALSE] / outer(s, s)
  R <- tryCatch(chol(C), error = function(e) NULL)
  if (!is.null(R)) {
    L[on, ] <- s * t(R)
    return(L)
  }
  e <- eigen((C + t(C)) / 2, symmetric = TRUE)
  k <- which(e$values > max(e$values) * 1e-12)
  L <- matrix(0, nrow(V), length(k))
  L[on, ] <- s * (e$vectors[, k, drop = FALSE] %*%
                    diag(sqrt(e$values[k]), length(k)))
  L
}

#' A factor of `V[pos, pos]` of a sparse joint covariance held in whitened
#' coordinates (jcs_joint_cov()): `M[pos, ci] L`, with `V_gamma[ci, ci] =
#' L L'` over the whitened positions `ci` the rows load. The natural block
#' has the design's condition number squared, and `a' V a` formed from it
#' lost 9.7e-4 of a fitted value's variance where the factor keeps every
#' digit (frm_lp_basis(), frm_curve(); review of lane sparsese, F1).
#'
#' @noRd
jcs_root_at <- function(S, pos) {
  Mi <- S[["M"]][pos, , drop = FALSE]
  ci <- which(diff(Mi@p) > 0L)
  L <- jcs_psd_factor(jcs_block(S, ci, ci))
  out <- as.matrix(Mi[, ci, drop = FALSE] %*% L)
  rownames(out) <- S[["names"]][pos]
  out
}

#' fd_joint_root()'s factor on a sparse joint covariance held in whitened
#' coordinates: positions `pos`, the first `nd` the drawn parameters, the
#' rest group effects. `V[pos, pos] = Vrest + K K'`, with a column of `K`
#' per direction of the drawn parameters (`M` of their whitened block's
#' factor `Ld`, and `V_bd Ld'^+` on the group effects) and `Vrest` zero
#' but on the group effects, where it is their covariance given the drawn
#' parameters. fit_fd_se() differences the output along each column of
#' `K`, two evaluations each, and keeps its chain and batch routes for the
#' group effects; a factor with a column per group effect would cost two
#' evaluations per level, as the dense route's does (62 s on a cumulative
#' fit at 4508 joint parameters, review of lane sparsese). NULL where a
#' drawn row loads a group effect's position.
#'
#' @noRd
jcs_fd_root <- function(S, pos, nd) {
  d <- seq_len(nd)
  b <- setdiff(seq_along(pos), d)
  Md <- S[["M"]][pos[d], , drop = FALSE]
  cd <- which(diff(Md@p) > 0L)
  if (any(pos[b] %in% cd)) return(NULL)
  cg <- c(cd, pos[b])
  Vg <- jcs_block(S, cg, cg)
  k <- seq_along(cd)
  Ld <- jcs_psd_factor(Vg[k, k, drop = FALSE])
  Kb <- if (length(b)) {
    t(qr.solve(Ld, Vg[k, -k, drop = FALSE]))
  } else matrix(0, 0L, ncol(Ld))
  Vrest <- matrix(0, length(pos), length(pos))
  if (length(b)) Vrest[b, b] <- Vg[-k, -k, drop = FALSE] - tcrossprod(Kb)
  list(K = rbind(as.matrix(Md[, cd, drop = FALSE] %*% Ld), Kb),
       Vrest = Vrest)
}
