# spde(), a Matern field over a finite-element mesh, registered with
# core through frmtmb_register_covstruct() (R/zzz.R). Everything here
# moved from frmtmb's R/covstruct.R unchanged in arithmetic; the fits it
# gives are the fits core gave (dev/spatial-fitcompare.R in the frmtmb
# repository).

#' Matern fields over a mesh: `spde()`
#'
#' `spde(fem, gr = node)` in a formula adds a gaussian Matern field
#' through the SPDE representation of Lindgren, Rue and Lindstrom
#' (2011): on a fixed mesh the field's precision is
#' `Q = tau^2 (kappa^4 M0 + 2 kappa^2 M1 + M2)`, with the mesh's three
#' finite-element matrices as fixed data. It is a predictor term, not a
#' bar term, and it is available once frmtmb.spatial is loaded:
#' `library(frmtmb.spatial)`. Without the package, [frmtmb::frm()]
#' refuses the term and names this package. spde() is frmtmb's own
#' spelling; brms has no such term.
#'
#' @section Arguments of the term:
#' \describe{
#'   \item{`fem`}{The mesh's finite-element matrices, as a list holding
#'     `c0`, `g1` and `g2` (`fmesher::fm_fem()`'s names) or `M0`, `M1`
#'     and `M2` (INLA's). They must be square and of one size, one row
#'     per mesh node. Put the list in `data2 = list(fem = fem)`. Mesh
#'     construction stays outside the package. `alpha = 2`, the order
#'     the three-matrix form covers, is the only one.}
#'   \item{`gr`}{The node of each observation, by ROW NUMBER: whole
#'     numbers in `1..nrow(M0)`, as an integer or numeric column, a
#'     character column of those numbers, or a factor whose levels are
#'     those numbers. The matrices carry no dimnames, so nothing else
#'     can reconcile a set of labels with them, and any other labeling
#'     is refused. A node no row visits keeps its column and its prior.}
#' }
#'
#' @section What the fit reports:
#' [frmtmb::confint_varcorr()] reports the marginal standard deviation
#' `sd(spde)` and the range `range(spde)`, under the planar
#' `alpha = 2` identities `range = sqrt(8) / kappa` and
#' `sigma^2 = 1 / (4 pi kappa^2 tau^2)`. Every mesh node is a level of
#' the block, so newdata names a node by its row number too. A general
#' projector matrix, which spreads an observation over the vertices of
#' the triangle it falls in, is not supported yet.
#'
#' @return `spde()` is a formula term, not a function: [frmtmb::frm()]
#'   reads it at parse time, and the value it contributes is a
#'   random-effect block of the model. This page documents the term
#'   and returns nothing.
#' @seealso `vignette("spatial", package = "frmtmb.spatial")`, and
#'   [frmtmb::frmtmb_register_covstruct()] for the seam the term is
#'   registered through.
#' @examples
#' # the alpha = 2 finite-element triple of a regular 1-D mesh
#' nn <- 10
#' h <- 0.5
#' C0 <- diag(rep(h, nn))
#' C0[1, 1] <- C0[nn, nn] <- h / 2
#' G <- matrix(0, nn, nn)
#' for (i in seq_len(nn - 1)) {
#'   G[i, i] <- G[i, i] + 1 / h
#'   G[i + 1, i + 1] <- G[i + 1, i + 1] + 1 / h
#'   G[i, i + 1] <- G[i + 1, i] <- -1 / h
#' }
#' fem <- list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
#' set.seed(5)
#' d <- data.frame(node = rep(seq_len(nn), each = 5), x = rnorm(5 * nn))
#' d$y <- 0.3 + 0.4 * d$x + sin(d$node / 2) + rnorm(nrow(d), 0, 0.3)
#' fit <- frm(bf(y ~ x + spde(fem, gr = node)), family = gaussian(),
#'            data = d, data2 = list(fem = fem))
#' confint_varcorr(fit)
#' @name frmtmb-spde
#' @aliases spde
NULL

# ------------------------------------------------------------------ SPDE
#
# Gaussian Matern field through the SPDE / finite-element
# representation (Lindgren-Rue-Lindstrom): on a fixed mesh the field's
# precision is
#
#   Q(kappa, tau) = tau^2 (kappa^4 M0 + 2 kappa^2 M1 + M2)
#
# with M0/M1/M2 (fmesher's c0/g1/g2, INLA's M0/M1/M2) the mesh's finite
# element matrices - fixed sparse data. Only the three weights depend on
# the parameters, so the whole precision is an AD-weighted sum of fixed
# sparse matrices and the log-determinant comes from RTMB's sparse GMRF
# density; theta = (log tau, log kappa). alpha = 2 (nu = 1 in the plane)
# is the only order the three-matrix form covers, which is the standard
# choice and what fm_fem() returns.

#' The three finite-element matrices, under either the INLA (M0/M1/M2)
#' or the fmesher (c0/g1/g2) spelling. The mesh size is whatever the
#' matrices say it is - see spde_node_index() for why the data cannot be
#' asked instead.
#'
#' @noRd
spde_matrices <- function(fem) {
  nms <- names(fem) %||% character(0)
  key <- if (all(c("M0", "M1", "M2") %in% nms)) c("M0", "M1", "M2")
         else if (all(c("c0", "g1", "g2") %in% nms)) c("c0", "g1", "g2")
         else NULL
  if (is.null(key)) {
    frm_stop("spde(): fem must be a list holding M0, M1, M2 (INLA) or ",
             "c0, g1, g2 (fmesher::fm_fem)", call. = FALSE)
  }
  M0 <- fem[[key[1L]]]
  if (length(dim(M0)) != 2L || nrow(M0) != ncol(M0)) {
    frm_stop("spde(): ", key[1L], " must be a square matrix, one row per ",
             "mesh node", call. = FALSE)
  }
  n <- nrow(M0)
  out <- list()
  for (i in seq_along(key)) {
    Mi <- fem[[key[i]]]
    if (length(dim(Mi)) != 2L || nrow(Mi) != n || ncol(Mi) != n) {
      frm_stop("spde(): ", key[i], " is ",
               paste(dim(Mi) %||% length(Mi), collapse = " x "), " but ",
               key[1L], " has ", n, " mesh nodes; the three finite-element ",
               "matrices must be square and the same size", call. = FALSE)
    }
    out[[c("M0", "M1", "M2")[i]]] <-
      methods::as(methods::as(Matrix::Matrix(Mi, sparse = TRUE),
                              "generalMatrix"), "CsparseMatrix")
  }
  out
}

#' Node indices for an spde() grouping variable.
#'
#' The finite-element matrices are indexed by mesh ROW NUMBER and carry
#' no dimnames to match labels against, so nothing reconciles a level
#' ordering with the mesh. Deriving the ordering from the grouping
#' variable - factor levels, or sort(unique(as.character(gv))) - reads
#' integer node ids lexicographically ("1", "10", "11", "2", ...), which
#' permutes the field against its own precision and converges silently to
#' a different model. The contract is therefore explicit: gr holds mesh
#' row indices, in 1..nrow(M0). Anything else is refused rather than
#' guessed at. Unobserved nodes are fine (they keep their prior); gaps
#' are not an error.
#'
#' @noRd
spde_node_index <- function(gv, n, gr_expr) {
  idx <- spde_node_numeric(gv)
  if (is.null(idx) || any(idx != trunc(idx))) {
    frm_stop("spde(): mesh nodes are indexed by their ROW NUMBER in the ",
             "finite-element matrices, which carry no location names, so ",
             "gr must hold whole-number node indices in 1..", n, "; '",
             deparse1(gr_expr), "' holds labels that are not whole numbers. ",
             "Map each observation onto its mesh row first and pass those ",
             "indices as gr", call. = FALSE)
  }
  idx <- as.integer(idx)
  bad <- unique(idx[idx < 1L | idx > n])
  if (length(bad)) {
    frm_stop("spde(): gr = ", deparse1(gr_expr), " holds node index/indices ",
             paste(utils::head(sort(bad), 5), collapse = ", "),
             ", outside the mesh's 1..", n, " rows", call. = FALSE)
  }
  idx
}

#' The numeric node ids behind an spde() grouping variable, or NULL when
#' the values are not numbers. A factor is read through its LEVELS, never
#' its integer codes: the codes are a level ordering, which is exactly
#' the thing that must not decide mesh rows.
#'
#' @noRd
spde_node_numeric <- function(gv) {
  v <- if (is.factor(gv)) {
    lv <- suppressWarnings(as.numeric(levels(gv)))
    if (anyNA(lv)) return(NULL)
    lv[as.integer(gv)]
  } else if (is.numeric(gv)) {
    as.numeric(gv)
  } else {
    w <- suppressWarnings(as.numeric(as.character(gv)))
    if (anyNA(w)) return(NULL)
    w
  }
  if (anyNA(v)) NULL else v
}

#' The block level labels an spde() grouping variable maps onto, used by
#' both frame assembly and newdata prediction so the two spellings of a
#' node ("3", 3L, factor "3") always land on the same column.
#'
#' @noRd
spde_node_labels <- function(gv) {
  v <- spde_node_numeric(gv)
  if (is.null(v) || any(v != trunc(v))) return(as.character(gv))
  as.character(as.integer(v))
}

#' Number of `theta` entries for an spde block: `log tau` and `log kappa`.
#'
#' @noRd
spde_npar <- function() 2L

#' Starting `theta` for an spde block: unit `tau` and unit `kappa`.
#'
#' @noRd
spde_start <- function() c(0, 0)

# Range and marginal sd follow the alpha = 2, d = 2 identities
# (Lindgren et al. 2011 eq. 2): range = sqrt(8 nu)/kappa with nu = 1,
# sigma^2 = 1 / (4 pi kappa^2 tau^2).
#' Marginal range of an spde field at a `theta`, on the natural scale.
#'
#' @noRd
spde_range <- function(theta) sqrt(8) / exp(theta[2])

#' Marginal standard deviation of an spde field at a `theta`, on the
#' natural scale.
#'
#' @noRd
spde_sd <- function(theta) {
  1 / (exp(theta[1]) * exp(theta[2]) * sqrt(4 * pi))
}

spde_slots <- list(
  npar = function(blk) spde_npar(),
  start = function(blk) spde_start(),
  lpdf = function(b, theta, blk) {
    a <- blk[["aux"]]
    kap2 <- exp(2 * theta[2])
    Q <- exp(2 * theta[1]) * (kap2 * kap2 * a[["M0"]] +
                                (2 * kap2) * a[["M1"]] + a[["M2"]])
    sum(RTMB::dgmrf(b, 0, Q, log = TRUE))
  },
  vcov = function(theta, blk) {
    matrix(spde_sd(theta)^2, 1, 1,
           dimnames = list("sd(spde)", "sd(spde)"))
  },
  # sigma and range are analytic functions of (log tau, log kappa):
  # log sigma = -log tau - log kappa - log(4 pi) / 2 and
  # log range = log(8) / 2 - log kappa, so the delta method is one
  # exact linear map with no differencing
  natural = function(theta, blk) {
    list(term = c("sd(spde)", "range(spde)"), type = c("sd", "range"),
         est = c(log(spde_sd(theta)), log(spde_range(theta))),
         grad = rbind(c(-1, -1), c(0, -1)))
  },
  draw = function(theta, blk) {
    U <- Matrix::chol(spde_prec(theta, blk))
    as.vector(Matrix::solve(U, stats::rnorm(blk[["n_levels"]])))
  },
  levels_coupled = TRUE,
  term = list(
    args = c("fem", "gr"),
    usage = "spde(fem, gr = g)",
    parse = function(args, env) {
      if (is.null(args[["fem"]]) || is.null(args[["gr"]])) {
        frm_stop("spde() needs the mesh matrices and a grouping variable: ",
                 "spde(fm_fem(mesh), gr = node)", call. = FALSE)
      }
      list(fem_expr = args[["fem"]], gr_expr = args[["gr"]])
    },
    build = function(spec, gr, lookup) spde_term_build(spec, gr, lookup),
    # the block's levels are mesh ROW NUMBERS, so a node is read as a
    # number on newdata too: as.character() on a double would spell
    # node 100000 as "1e+05" and lose the column
    levels_of = function(gr) spde_node_labels(gr),
    # every mesh node is a level, observed or not, so a node outside
    # the fitted set is not an unseen level of anything
    levels_complete = TRUE
  )
)

#' The spde() block over a mesh: the mesh names nothing, so the block's
#' levels ARE the mesh rows and the data is permuted to them. Every node
#' gets a column, observed or not (an unobserved one keeps its prior).
#'
#' @noRd
spde_term_build <- function(spec, gr, lookup) {
  aux <- spde_matrices(lookup(spec[["fem_expr"]], "spde()"))
  n_node <- nrow(aux[["M0"]])
  list(levels = as.character(seq_len(n_node)),
       index = spde_node_index(gr, n_node, spec[["gr_expr"]]),
       aux = aux)
}

#' Numeric precision of an spde block at a theta (draws, diagnostics).
#'
#' @noRd
spde_prec <- function(theta, blk) {
  a <- blk[["aux"]]
  kap2 <- exp(2 * theta[2])
  exp(2 * theta[1]) * (kap2 * kap2 * a[["M0"]] + (2 * kap2) * a[["M1"]] +
                         a[["M2"]])
}


#' Register spde() with frmtmb. Called from .onLoad().
#'
#' @noRd
register_spde <- function() {
  do.call(frmtmb_register_covstruct, c(list(name = "spde"), spde_slots))
}
