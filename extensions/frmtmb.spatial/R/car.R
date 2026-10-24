# car(), brms's conditional autoregressive field, registered with core
# through frmtmb_register_covstruct() (R/zzz.R). Everything here moved
# from frmtmb's R/covstruct.R, R/parse.R and R/frame.R unchanged in
# arithmetic; the fits it gives are the fits core gave
# (dev/spatial-fitcompare.R in the frmtmb repository).

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Conditional autoregressive fields: `car()`
#'
#' `car(M, gr = g, type = )` in a formula adds a conditional
#' autoregressive (CAR) field over the locations an adjacency matrix
#' names, with brms's spelling and its four types. It is a predictor
#' term, not a bar term, and it is available once frmtmb.spatial is
#' loaded: `library(frmtmb.spatial)`. Without the package, [frmtmb::frm()]
#' refuses the term and names this package.
#'
#' @section Arguments of the term:
#' \describe{
#'   \item{`M`}{A symmetric adjacency matrix whose dimnames (row names,
#'     column names, or both, which then have to agree) name the
#'     locations, with a row for every location the data hold. Entries
#'     must be present and non-negative; non-zero weights are read as 1,
#'     with a message. Put it in `data2 = list(M = M)`.}
#'   \item{`gr`}{The grouping variable that maps observations onto
#'     locations. It must be named: brms's deprecated `gr = NA` is
#'     refused.}
#'   \item{`type`}{`"escar"` (the default) is the proper CAR,
#'     `Q = tau (D - rho W)`, and needs every location to have a
#'     neighbor. `"esicar"` is the intrinsic CAR under an EXACT
#'     sum-to-zero constraint, `"icar"` the same field under brms's SOFT
#'     constraint at scale `con_sd`, and `"bym2"` the scaled mixture of a
#'     spatial and an unstructured part.}
#'   \item{`con_sd`}{The scale of the soft constraint, brms's `1e-3` by
#'     default. It changes nothing an `"esicar"` term reports.}
#' }
#'
#' @section What the fit reports:
#' [frmtmb::confint_varcorr()] reports `sd(car)`, the field scale, and
#' for `"escar"` and `"bym2"` the mixing parameter on `(0, 1)` under
#' brms's names `car` and `rhocar`. [frmtmb::ranef()] gives the field at
#' each location; for `"esicar"` that is the centered field, with its
#' conditional standard deviations. A location the fit did not see has
#' no law of its own, so `predict()` with `allow_new_levels =
#' TRUE` predicts it at the population level.
#'
#' Every constrained type constrains each connected component of the
#' graph, where brms constrains the global sum; on a connected graph the
#' two are the same constraint.
#'
#' @return `car()` is a formula term, not a function: [frmtmb::frm()]
#'   reads it at parse time, and the value it contributes is a
#'   random-effect block of the model. This page documents the term
#'   and returns nothing.
#' @seealso `vignette("spatial", package = "frmtmb.spatial")`, and
#'   [frmtmb::frmtmb_register_covstruct()] for the seam the term is
#'   registered through.
#' @examples
#' set.seed(4)
#' grid <- expand.grid(r = 1:4, c = 1:4)
#' n <- nrow(grid)
#' W <- outer(seq_len(n), seq_len(n), function(i, j) {
#'   as.numeric(abs(grid$r[i] - grid$r[j]) + abs(grid$c[i] - grid$c[j]) == 1)
#' })
#' dimnames(W) <- list(paste0("L", 1:n), paste0("L", 1:n))
#' d <- data.frame(loc = factor(rep(rownames(W), each = 6),
#'                              levels = rownames(W)))
#' d$y <- 1 + rnorm(n)[as.integer(d$loc)] + rnorm(nrow(d), 0, 0.5)
#' fit <- frm(bf(y ~ car(W, gr = loc, type = "icar")), family = gaussian(),
#'            data = d, data2 = list(W = W))
#' confint_varcorr(fit)
#' @name frmtmb-car
#' @aliases car
NULL

# ------------------------------------------------------- CAR / ICAR / BYM2
#
# Spatial conditional-autoregressive fields over a fixed adjacency
# matrix - brms's car(M, gr = g, type = ). The field lives on the
# grouping levels (one intercept per location), so the block is
# dim = 1 with n_levels locations, exactly like gr(prec = Q); what
# car() adds is that the precision is BUILT from the graph and carries
# its own hyperparameters.
#
# Parameterization follows brms: theta = log(sdcar) for the intrinsic
# forms and (log sdcar, logit rho) for the two that mix, with
# tau = 1 / sdcar^2 the precision multiplier.
#
#   escar    Q = tau (D - rho W), proper for rho in (0, 1); no
#            constraint, because a proper CAR needs none
#   esicar   Q = tau (D - W) on the sum-to-zero subspace, the
#            constraint imposed EXACTLY
#   icar     the same intrinsic precision under a SOFT sum-to-zero
#            constraint, brms's normal(sum(zcar) | 0, con_sd Nloc)
#   bym2     sd^2 [(1 - rho) I + (rho / scale) K^-1], the
#            Riebler et al. scaled mixture brms implements
#
# esicar and icar are two different models, as they are in brms. The
# two spellings selected one density here through 0.51.0; what each one
# now costs is below.
#
# THE NORMALIZING CONSTANT is analytic in every case, so no
# normalize-trick / on-tape log-determinant is needed and the density
# costs one sparse matrix-vector product:
#
#   escar   log|Q| = n log tau + sum_i log d_i + sum_i log(1 - rho e_i)
#           with e_i the eigenvalues of D^-1/2 W D^-1/2 (fixed data).
#   icar    the graph Laplacian L = D - W is rank n - c with c the
#           number of connected components, so log|Q|* = (n - c) log
#           tau + log|L|*. The soft constraint below turns that into an
#           exact n log tau.
#   esicar  the pseudo-determinant itself, (n - c) log tau + log|L|*,
#           which is what brms's sparse_icar_lpdf computes.
#
# THE SUM-TO-ZERO CONSTRAINT. An intrinsic CAR is improper: L
# annihilates the indicator s_j of each connected component, so shifting
# the field inside a component and the intercept the other way leaves
# the likelihood untouched and the ML problem is rank deficient, not
# merely ill conditioned. The two intrinsic types answer that
# differently, and both answers are brms's.
#
# icar, the SOFT answer. A rank-c update rides on tau:
#
#   Q(tau) = tau * K,   K = L + sum_j kappa_j s_j s_j',
#   kappa_j = 1 / (con_sd n_j)^2,
#
# with n_j the size of component j. K is fixed data, so
# log|Q| = n log tau + log|K| is exact and constant-free. The component
# sums are pinned at an sd of con_sd n_j sdcar rather than exactly
# zero, which is the same model as an extra random intercept of sd
# con_sd sdcar per component. That intercept is confounded with beta_0,
# and its convolution IS the bias against the exact constraint:
# measured on a 4 x 4 lattice, 1e-3 is off by 4.7e-4 in the
# log-likelihood (3.6e-5 relative in sdcar), 1e-4 by 4.7e-6 (3.6e-7),
# 1e-5 by 4.5e-8 (3.7e-9), 1e-6 by 9.2e-10, and 1e-7 loses to roundoff
# (1.1e-4). Tightening also costs optimizer robustness. The constraint
# direction carries a factor con_sd^-2 of the block Hessian, and over 25
# lattice refits nlminb reported false convergence 0 times at 1e-3, once
# at 1e-4 and 6 times at 1e-5, so the loose default is the better trade.
# con_sd defaults to brms's 1e-3, so the same call is the same model
# here and there.
#
# esicar, the EXACT answer, and it never solves a constrained problem.
# The field that reaches the linear predictor is the CENTERED b,
#
#   c = b - sum_j s_j (s_j'b / n_j),
#
# applied in expand_b(), so 1'c = 0 within every component by
# construction. The component means m_j = s_j'b / n_j are then inert:
# they enter no linear predictor and no cross term. Give them a
# tau-FREE density and the block precision is
#
#   Q(tau) = tau L + P0,   P0 = sum_j kappa_j s_j s_j',
#
# which is block diagonal on the orthogonal split R^n = range(L) plus
# span{s_j}, because L kills every s_j and P0 kills everything
# orthogonal to them. Three consequences, and they are the design:
#
#   1. log|Q| = (n - c) log tau + log|K| with the SAME K icar builds,
#      since |tau L + P0| and |L + P0| differ only on range(L). The
#      pseudo-determinant is exact and costs no new factorization.
#   2. p(b) = p_hard(c) * prod_j N(m_j; 0, con_sd^2) exactly, and the
#      predictor depends on c alone, so the m integral is 1 and the
#      marginal likelihood IS the hard-constrained one. The Laplace
#      approximation splits the same way, with an exactly Gaussian m
#      block, so this holds to machine precision and not in a limit.
#   3. logLik() is INVARIANT to con_sd, which is the signature that
#      separates the two types: icar's moves with it and esicar's does
#      not. con_sd only scales an inert coordinate here, and it is left
#      at kappa_j so that the block Hessian's conditioning is icar's.
#
# m* = 0 exactly at the mode, because the m block is quadratic and
# decoupled and one Newton step lands on it. So `b` is already centered
# wherever anything reads it, and the inert coordinate carries variance
# con_sd^2 = 1e-6, four orders below sdcar^2 on any real fit. That is
# why the conditional-mode standard errors need no branch of their own.
#
# brms constrains the GLOBAL sum for esicar (rcar[Nloc] = -sum(zcar))
# and normalizes by (Nloc - 1) log tau, which is right for a connected
# graph and improper otherwise. The per-component form above agrees
# with brms whenever c = 1 and stays proper when it is not.
#
# The price of any sum-to-zero constraint is a dense rank-c update
# inside the block's Laplace Hessian, which caps the practical field
# size in the low thousands.
car_con_sd_default <- 1e-3
car_types <- c("escar", "esicar", "icar", "bym2")

#' brms's validate_car_matrix plus the level matching data_ac() does: a
#' symmetric binary adjacency over exactly the locations the data show,
#' in the block's level order.
#'
#' @noRd
car_adjacency <- function(M, locs) {
  if (length(dim(M)) != 2L || nrow(M) != ncol(M)) {
    frm_stop("car(): M must be a square adjacency matrix", call. = FALSE)
  }
  Md <- as.matrix(M)
  # NA first: every check below compares, and a comparison against NA is
  # NA, which reaches the user as "missing value where TRUE/FALSE
  # needed" from isSymmetric() - naming neither the matrix nor the cell
  if (anyNA(Md)) {
    bad <- which(is.na(Md), arr.ind = TRUE)
    frm_stop("car(): M has ", nrow(bad), " missing entry/entries (the first ",
             "at row ", bad[1L, 1L], ", column ", bad[1L, 2L], "); an ",
             "adjacency matrix needs a 0 for every non-neighbor pair",
             call. = FALSE)
  }
  # Locations are matched by NAME, so one set of names is enough - but
  # when both are present they have to agree, or the row a location gets
  # and the column it gets are two different places.
  rn <- rownames(Md)
  cn <- colnames(Md)
  if (is.null(rn) && is.null(cn)) {
    frm_stop("car(): M needs dimnames naming the locations", call. = FALSE)
  }
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
    frm_stop("car(): M's rownames and colnames must name the same locations ",
             "in the same order", call. = FALSE)
  }
  nms <- rn %||% cn
  dup <- unique(nms[duplicated(nms)])
  if (length(dup)) {
    frm_stop("car(): M names location(s) ",
             paste(utils::head(dup, 5), collapse = ", "),
             " more than once", call. = FALSE)
  }
  miss <- setdiff(locs, nms)
  if (length(miss)) {
    frm_stop("car(): M has no row for location(s) ",
             paste(utils::head(miss, 5), collapse = ", "),
             if (length(miss) > 5) " ..." else "", call. = FALSE)
  }
  # positional subsetting, so a matrix carrying only rownames (or only
  # colnames) is read exactly like one carrying both
  pos <- match(locs, nms)
  W <- methods::as(Matrix::Matrix(Md[pos, pos, drop = FALSE],
                                  sparse = TRUE), "generalMatrix")
  if (!Matrix::isSymmetric(W, check.attributes = FALSE)) {
    frm_stop("car(): M must be symmetric", call. = FALSE)
  }
  if (any(Matrix::diag(W) != 0)) {
    frm_stop("car(): M must have a zero diagonal (no location neighbors ",
             "itself)", call. = FALSE)
  }
  W <- methods::as(W, "CsparseMatrix")
  # brms's validate_car_matrix reads an adjacency as non-negative
  # weights. Binarizing a negative entry would turn a stated repulsion
  # into a neighbor, so refuse it rather than reinterpret it.
  if (any(W@x < 0)) {
    frm_stop("car(): M has ", sum(W@x < 0), " negative entry/entries (the ",
             "smallest is ", format(min(W@x), digits = 4), "); an adjacency ",
             "matrix holds non-negative weights", call. = FALSE)
  }
  if (any(W@x != 1)) {
    frm_message("car(): converting all non-zero values in M to 1.")
    W@x[W@x != 1] <- 1
  }
  dimnames(W) <- list(locs, locs)
  W
}

#' Number of `theta` entries a CAR block needs. The two types that mix an
#' independent part with the spatial one carry a `rho`; the intrinsic
#' forms do not.
#'
#' @noRd
car_npar <- function(type) if (type %in% c("escar", "bym2")) 2L else 1L

#' Starting `theta` for a CAR block: unit `sdcar`, and `rho` at the middle
#' of its range where the type has one.
#'
#' @noRd
car_start <- function(type) numeric(car_npar(type))

# brms bounds both mixing parameters on (0, 1)
#' Maps an unbounded `theta` entry onto a CAR mixing parameter in (0, 1).
#'
#' @noRd
car_rho <- function(x) 1 / (1 + exp(-x))

#' Connected components of a symmetric adjacency matrix, by breadth-first
#' sweep over the sparse column pattern (isolated levels are their own
#' component, which is what the rank correction counts).
#'
#' @noRd
car_components <- function(W) {
  n <- nrow(W)
  Wt <- methods::as(W, "TsparseMatrix")
  nbr <- split(c(Wt@j + 1L, Wt@i + 1L),
               factor(c(Wt@i + 1L, Wt@j + 1L), levels = seq_len(n)))
  comp <- integer(n)
  k <- 0L
  for (s in seq_len(n)) {
    if (comp[s]) next
    k <- k + 1L
    queue <- s
    comp[s] <- k
    while (length(queue)) {
      v <- queue[1L]
      queue <- queue[-1L]
      nb <- nbr[[v]]
      new <- nb[comp[nb] == 0L]
      if (length(new)) {
        comp[new] <- k
        queue <- c(queue, unique(new))
      }
    }
  }
  comp
}

#' brms's scaling factor for BYM2 (brms:::.car_scale): the geometric mean
#' of the marginal variances of a unit-scale ICAR field under the
#' sum-to-zero constraint. Reproduced here (perturbation included) so the
#' rho of a bym2 fit means what it means in brms.
#'
#' @noRd
car_scale_factor <- function(W) {
  n <- nrow(W)
  Q <- Matrix::Diagonal(n, Matrix::rowSums(W)) - W
  Q <- Q + Matrix::Diagonal(n) * max(Matrix::diag(Q)) *
    sqrt(.Machine$double.eps)
  Sigma <- Matrix::solve(Q)
  A <- matrix(1, 1, n)
  Wc <- Sigma %*% t(A)
  Sigma <- Sigma - Wc %*% solve(A %*% Wc) %*% Matrix::t(Wc)
  exp(mean(log(Matrix::diag(Sigma))))
}

#' Everything the density needs, precomputed once from the adjacency
#' matrix over the grouping levels. `type` decides which fields matter;
#' all of them are fixed data.
#'
#' @noRd
car_aux <- function(W, type, con_sd = car_con_sd_default) {
  n <- nrow(W)
  deg <- as.numeric(Matrix::rowSums(W))
  # Csparse, not just generalMatrix: `L %*% b` runs on the tape for
  # esicar, and RTMB's sparse product wants the compressed form
  L <- methods::as(methods::as(Matrix::Diagonal(n, deg) - W,
                               "generalMatrix"), "CsparseMatrix")
  aux <- list(type = type, n = n, W = W, L = L, deg = deg)
  if (type == "escar") {
    if (any(deg == 0)) {
      frm_stop("car(type = \"escar\"): every location needs at least one ",
               "neighbor; location(s) ",
               paste(rownames(W)[deg == 0], collapse = ", "),
               " have none. Use type = \"icar\" instead", call. = FALSE)
    }
    isq <- diag(1 / sqrt(deg), nrow = n)
    ev <- eigen(isq %*% as.matrix(W) %*% isq, symmetric = TRUE,
                only.values = TRUE)$values
    # The normalized adjacency's spectrum lies in [-1, 1], but LAPACK
    # returns the extreme eigenvalues a few ulp outside it. The density
    # takes log(1 - rho * eig) with rho = plogis(theta2), which reaches
    # 1 to double precision by theta2 = 36: a stray eig = 1 + 2e-16
    # makes the argument negative and the whole objective NaN, and the
    # optimizer only reports "NA/NaN function evaluation". Clamping
    # costs nothing elsewhere - the determinant is unchanged to twelve
    # digits - and keeps the boundary evaluable.
    aux$eigW <- pmin(pmax(ev, -1), 1)
    aux$ldet_deg <- sum(log(deg))
    return(aux)
  }
  comp <- car_components(W)
  nj <- as.numeric(table(comp))
  aux$con_sd <- con_sd
  aux$kappa0 <- 1 / (con_sd * nj)^2
  aux$Sgrp <- Matrix::sparseMatrix(i = comp, j = seq_len(n), x = 1,
                                   dims = c(length(nj), n))
  K <- L + Matrix::t(aux$Sgrp) %*% Matrix::Diagonal(length(nj),
                                                    aux$kappa0) %*%
    aux$Sgrp
  aux$K <- methods::as(methods::as(K, "generalMatrix"), "CsparseMatrix")
  aux$ldet_K <- as.numeric(Matrix::determinant(aux$K,
                                               logarithm = TRUE)$modulus)
  aux$n_comp <- length(nj)
  aux$comp <- comp
  aux$nj <- nj
  if (type == "bym2") {
    # the scaled mixture needs the covariance of the unit-scale ICAR
    aux$Kinv <- as.matrix(Matrix::solve(aux$K))
    aux$scale <- car_scale_factor(W)
  }
  aux
}

#' The esicar field: `b` with each connected component's mean removed,
#' so that the vector the linear predictor sees sums to zero within
#' every component by construction. Runs on the tape (advector `b`) and
#' off it, because expand_b() is the numeric path too.
#'
#' @noRd
car_center <- function(b, a) {
  m <- as.vector(a[["Sgrp"]] %*% b) / a[["nj"]]
  b - m[a[["comp"]]]
}

#' The derivative of `car_center()`. `d(field)/d(b)` is the
#' per-component centering projection
#' `P = I - sum_j (1 / n_j) s_j s_j'`, returned as triplets in
#' WITHIN-BLOCK positions for the caller to map through `c_idx` and
#' `b_idx`.
#'
#' It lives beside `car_center()` because the two have to agree. While
#' the delta method paired the `Z` columns with `b` through the
#' identity, every prediction standard error and every `ranef()`
#' conditional SD carried the component means - a coordinate the
#' predictor never sees and the data never speak to - as exactly
#' `con_sd^2` in the variance.
#'
#' A singleton component contributes NO entries: `b_i - b_i` is
#' identically zero, so its row and column of `P` are structurally
#' empty and its field is a constant, not a random variable.
#'
#' @noRd
car_center_jacobian <- function(a) {
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (mem in split(seq_len(a[["n"]]), a[["comp"]])) {
    k <- length(mem)
    if (k == 1L) next
    ii <- c(ii, rep(mem, times = k))
    jj <- c(jj, rep(mem, each = k))
    xx <- c(xx, as.vector(diag(1, k) - 1 / k))
  }
  list(i = ii, j = jj, x = xx)
}

#' Conditional SDs of the esicar FIELD, from the variances of `b` that
#' `sdreport()` returns. `ranef()` displays `P b`, so it has to report
#' `diag(P V P')` rather than `diag(V)`.
#'
#' In closed form rather than by materializing `V`. The block precision
#' `tau L + P0` is block diagonal on the orthogonal split
#' `b = f + sum_j m_j s_j`, and `m` enters neither the data term nor any
#' cross term, so `V = cov(f) + sum_j con_sd^2 s_j s_j'` with the two
#' terms orthogonal, `P` annihilates the second exactly, and
#'
#'     diag(P V P') = diag(V) - con_sd^2
#'
#' at every level of every component. `n_j` does not enter, because
#' `var(m_j)` is `con_sd^2` whatever the component's size, so a
#' singleton lands on exactly zero, which is the variance its
#' identically-zero field must have. `pmax()` is there for that case
#' only: the subtraction is exact to rounding, and the singleton's
#' difference can land a few ulp below zero.
#'
#' @noRd
car_center_condsd <- function(v, a) {
  sqrt(pmax(v - a[["con_sd"]]^2, 0))
}

#' The esicar block precision at a theta, `tau L + P0`. Unlike icar's
#' `tau K` the rank-c term does NOT ride on tau: that is what makes the
#' inert component means tau-free and the marginal likelihood the
#' hard-constrained one.
#'
#' @noRd
car_esicar_prec <- function(theta, a) {
  P0 <- Matrix::t(a[["Sgrp"]]) %*%
    Matrix::Diagonal(a[["n_comp"]], a[["kappa0"]]) %*% a[["Sgrp"]]
  exp(-2 * theta[1]) * a[["L"]] + P0
}

#' Numeric covariance of the whole field at a theta (draws, VarCorr
#' details); off the tape, so a dense solve is fine.
#'
#' @noRd
car_cov <- function(theta, blk) {
  a <- blk[["aux"]]
  s2 <- exp(2 * theta[1])
  if (a[["type"]] == "escar") {
    rho <- car_rho(theta[2])
    Q <- Matrix::Diagonal(a[["n"]], a[["deg"]]) - rho * a[["W"]]
    return(s2 * as.matrix(Matrix::solve(Q)))
  }
  if (a[["type"]] == "bym2") {
    rho <- car_rho(theta[2])
    return(s2 * ((1 - rho) * diag(a[["n"]]) +
                   (rho / a[["scale"]]) * a[["Kinv"]]))
  }
  if (a[["type"]] == "esicar") {
    # the covariance of `b`, which is what a draw needs; expand_b()
    # centers it afterwards and the centered draw is the constrained
    # field. Not s2 * anything: the rank-c term is tau-free.
    return(as.matrix(Matrix::solve(car_esicar_prec(theta, a))))
  }
  s2 * as.matrix(Matrix::solve(a[["K"]]))
}

car_slots <- list(
  sd_idx = 1L,
  lpdf = function(b, theta, blk) {
    a <- blk[["aux"]]
    n <- a[["n"]]
    if (a[["type"]] == "bym2") {
      # the mixture has no sparse precision, so the dense marginal
      # covariance is the honest form; one MVN observation
      rho <- car_rho(theta[2])
      Sigma <- RTMB::matrix(
        (1 - rho) * exp(2 * theta[1]) * as.vector(diag(n)) +
          (rho * exp(2 * theta[1]) / a[["scale"]]) *
          as.vector(a[["Kinv"]]), n, n)
      return(sum(RTMB::dmvnorm(b, 0, Sigma, log = TRUE)))
    }
    tau <- exp(-2 * theta[1])
    if (a[["type"]] == "esicar") {
      # Q = tau L + P0. The two terms live on complementary invariant
      # subspaces, so the determinant splits into the pseudo-determinant
      # brms normalizes by, (n - c) log tau + log|L|*, plus the inert
      # means' own constant, whose sum with log|L|* IS log|K|.
      quad_l <- sum(b * as.vector(a[["L"]] %*% b))
      quad_p <- sum(a[["kappa0"]] * as.vector(a[["Sgrp"]] %*% b)^2)
      ldet <- (n - a[["n_comp"]]) * log(tau) + a[["ldet_K"]]
      return(0.5 * (ldet - n * log(2 * pi)) -
               0.5 * (tau * quad_l + quad_p))
    }
    if (a[["type"]] == "escar") {
      rho <- car_rho(theta[2])
      quad <- sum(a[["deg"]] * b^2) -
        rho * sum(b * as.vector(a[["W"]] %*% b))
      ldet <- n * log(tau) + a[["ldet_deg"]] +
        sum(log(1 - rho * a[["eigW"]]))
    } else {
      quad <- sum(b * as.vector(a[["K"]] %*% b))
      ldet <- n * log(tau) + a[["ldet_K"]]
    }
    0.5 * (ldet - n * log(2 * pi)) - 0.5 * tau * quad
  },
  vcov = function(theta, blk) {
    # sdcar, brms's reported scale; the field covariance is car_cov()
    matrix(exp(theta[1])^2, 1, 1,
           dimnames = list("sd(car)", "sd(car)"))
  },
  natural = function(theta, blk) {
    if (length(theta) < 2L) {
      return(list(term = "sd(car)", type = "sd", est = theta[1L],
                  grad = matrix(1, 1L, 1L)))
    }
    # brms's names for the two mixing parameters, both on (0, 1)
    nm <- if (identical(blk[["aux"]][["type"]], "bym2")) "rhocar" else "car"
    list(term = c("sd(car)", nm), type = c("sd", "prop"),
         est = theta[1:2], grad = diag(2L))
  },
  # the whole field is one draw from its (dense) covariance
  draw = function(theta, blk) {
    K <- car_cov(theta, blk)
    drop(crossprod(chol(K), stats::rnorm(nrow(K))))
  },
  levels_coupled = TRUE,
  expand = list(
    apply = function(b, blk) car_center(b, blk[["aux"]]),
    jacobian = function(blk) car_center_jacobian(blk[["aux"]]),
    condsd = function(v, blk) car_center_condsd(v, blk[["aux"]])
  ),
  term = list(
    args = c("M", "gr", "type", "con_sd"),
    usage = "car(M, gr = g)",
    parse = function(args, env) car_term_parse(args, env),
    build = function(spec, gr, lookup) car_term_build(spec, gr, lookup)
  ),
  npar = function(blk) car_npar(blk[["aux"]][["type"]]),
  start = function(blk) car_start(blk[["aux"]][["type"]])
)

#' brms's car(M, gr, type): M is an adjacency matrix (looked up in the
#' data or the formula environment, as gr(cov = A) is), gr the grouping
#' variable naming the locations, type one of brms's four. `args` holds
#' the call's arguments, unevaluated, matched by the parser.
#'
#' @noRd
car_term_parse <- function(args, env) {
  if (is.null(args[["M"]])) {
    frm_stop("car() needs an adjacency matrix and a grouping variable: ",
             "car(M, gr = g, type = \"icar\")", call. = FALSE)
  }
  # brms's gr = NA default (one location per observation) is deprecated
  # there ("Using CAR terms without a grouping factor is deprecated") and
  # never supported here: the field's levels come from the factor, so the
  # call has to name one. car(M) is that default, written by omission
  if (is.null(args[["gr"]]) || identical(args[["gr"]], NA) ||
        identical(args[["gr"]], quote(NA))) {
    frm_stop("car(): gr must name a grouping variable. brms's gr = NA ",
             "(one location per observation) is deprecated; build the ",
             "location factor and pass it", call. = FALSE)
  }
  type <- if (is.null(args[["type"]])) "escar" else {
    tp <- eval(args[["type"]], env)
    if (!is.character(tp) || length(tp) != 1L || !(tp %in% car_types)) {
      frm_stop("car(): type must be one of ",
               paste0("\"", car_types, "\"", collapse = ", "), call. = FALSE)
    }
    tp
  }
  con_sd <- if (is.null(args[["con_sd"]])) car_con_sd_default else {
    car_con_sd_arg(args[["con_sd"]], env)
  }
  list(M_expr = args[["M"]], type = type, con_sd = con_sd)
}

#' The value of car(con_sd = ), with the refusals core's own specials
#' give a malformed tuning value, in the same words.
#'
#' @noRd
car_con_sd_arg <- function(expr, env) {
  val <- tryCatch(eval(expr, env), error = function(e) {
    frm_stop("car(): cannot evaluate con_sd = ", deparse1(expr), ": ",
             conditionMessage(e), call. = FALSE)
  })
  if (length(val) != 1L) {
    frm_stop("car(): con_sd = ", deparse1(expr),
             " must be a single value (got length ", length(val), ")",
             call. = FALSE)
  }
  if (!is.numeric(val) || !all(is.finite(val))) {
    frm_stop("car(): con_sd = ", deparse1(expr),
             " must be finite and numeric", call. = FALSE)
  }
  if (any(val <= 0)) {
    frm_stop("car(): con_sd = ", deparse1(expr), " must be positive",
             call. = FALSE)
  }
  as.numeric(val)
}

#' The car() block over the data's locations: the adjacency matrix names
#' its locations, so the block's level order is the data's and the
#' matrix is permuted to it.
#'
#' @noRd
car_term_build <- function(spec, gr, lookup) {
  locs <- if (is.factor(gr)) levels(droplevels(gr)) else {
    sort(unique(as.character(gr)))
  }
  W <- car_adjacency(lookup(spec[["M_expr"]], "car()"), locs)
  list(levels = locs, index = match(as.character(gr), locs),
       aux = car_aux(W, spec[["type"]], spec[["con_sd"]]),
       # esicar's coefficients are its parameters CENTERED, so the
       # frame has to say that expand_b() must run
       expands = identical(spec[["type"]], "esicar"))
}


#' Register car() with frmtmb. Called from .onLoad().
#'
#' @noRd
register_car <- function() {
  do.call(frmtmb_register_covstruct, c(list(name = "car"), car_slots))
}
