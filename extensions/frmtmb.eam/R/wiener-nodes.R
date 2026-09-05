# Fixed quadrature nodes for the across-trial variability integrals.
#
# The nodes are built ONCE, when the family object is constructed, and
# they are doubles from then on. Nothing about them depends on a
# parameter: the node count is an argument of wiener(), the abscissas
# are constants of the interval (0, 1), and a parameter only ever
# rescales the interval they are mapped onto. That is what makes the
# integral tapeable at all - a node count or a node position that moved
# with a parameter would be a branch on a parameter, which an RTMB tape
# cannot record.

#' Gauss-Legendre nodes and weights on the unit interval.
#'
#' Golub and Welsch: the nodes of the n-point rule are the eigenvalues
#' of the symmetric tridiagonal Jacobi matrix of the Legendre
#' recurrence, and the weights are the squared first components of its
#' eigenvectors. Written out rather than taken from a package so that
#' this family adds no dependency for ten lines of linear algebra, and
#' so that the rule is reproducible without one.
#'
#' Gauss-Legendre rather than the midpoint or trapezoid rules the
#' uniform's flat weight might suggest: the integrand is not the uniform
#' density, it is the Wiener density seen through it, which is analytic
#' on the open interval. An n-point Gauss rule is exact for polynomials
#' of degree 2n-1 and converges geometrically on an analytic integrand,
#' where a composite midpoint rule converges at order n^-2. The
#' measured cost of that choice is in tests/testthat/test-variability.R.
#'
#' @return A list with `x` in (0, 1) and `w` summing to 1, so that
#'   `sum(w * f(x))` is the MEAN of f over the interval rather than its
#'   integral. The mean is what a uniform density asks for, and it makes
#'   the one-node case a plain evaluation at the midpoint.
#'
#' @noRd
ddm_gauss_legendre <- function(n) {
  n <- as.integer(n)
  if (n == 1L) return(list(x = 0.5, w = 1))
  k <- seq_len(n - 1L)
  b <- k / sqrt(4 * k * k - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- b
  J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  z <- e$values[ord]                  # nodes on (-1, 1)
  v <- e$vectors[1L, ord]
  list(x = 0.5 * (z + 1), w = v * v)  # weights already sum to 1
}

#' The node set a variability family carries.
#'
#' An integral that is switched off gets the ONE-node rule rather than
#' no rule, and the family passes it a width of zero. The density then
#' has one code path instead of four, and the switched-off case is not
#' an approximation of anything: a single node at the midpoint of a
#' zero-width interval is an evaluation at the midpoint, exactly.
#'
#' The drift never appears here. Its integral is Gaussian against an
#' exponential-quadratic and has a closed form, so it costs one square
#' root rather than a node loop, and there is no count to choose.
#'
#' @noRd
ddm_nodes <- function(variability, nodes) {
  rule <- function(nm) {
    if (nm %in% variability) {
      ddm_gauss_legendre(nodes[[nm]])
    } else {
      list(x = 0.5, w = 1)
    }
  }
  list(sz = rule("sz"), st = rule("st"))
}

#' Gauss-Hermite nodes and weights for a standard normal.
#'
#' The same Golub and Welsch construction on the probabilists' Hermite
#' recurrence, so `sum(w * f(x))` is the MEAN of f under `N(0, 1)`.
#'
#' Used only off the tape, by the post-fit mean. The density itself
#' never needs it: there the drift integral has a closed form, and a
#' node rule would be a worse answer to a solved problem. The mean
#' response time conditional on a boundary is a ratio of two integrals
#' over the drift and has no such form, which is why the nodes exist at
#' all.
#'
#' @noRd
ddm_gauss_hermite <- function(n) {
  n <- as.integer(n)
  if (n == 1L) return(list(x = 0, w = 1))
  k <- seq_len(n - 1L)
  b <- sqrt(k)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- b
  J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  v <- e$vectors[1L, ord]
  list(x = e$values[ord], w = v * v)
}
