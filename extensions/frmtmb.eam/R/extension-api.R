# The first of the two things this package offers ANOTHER extension
# package. R/ndt-seam.R carries the second, the BOUND a non-decision
# time is measured against, which item 1.0b asked for and which this
# file's closing paragraph anticipated.
#
# WHY IT EXISTS. frmtmb.learn's `rlddm()` is a delta learning rule whose
# value difference drives the drift rate of a Wiener first-passage
# density. The learning rule belongs there and the density belongs here,
# and until this file there was no exported route to the density: every
# Wiener function in this package is `@noRd`, so the only way to reach
# one was `frmtmb.eam:::ddm_lpdf_both()`, which is a promise nobody
# made. This is that promise, made deliberately and kept small.
#
# WHAT IS AND IS NOT PROMISED. The log density of the two-boundary
# Wiener first-passage time, with the parameterization `wiener()` uses
# and with the tape-safety `wiener()` needs. Nothing else here is
# public: the series truncations, the blend between them, the
# across-trial variability integrals and the CDF stay internal, and a
# caller who wants any of those should ask for a second seam rather than
# reach for a colon.

#' The Wiener first-passage log density, for another package's family
#'
#' The log density of the two-boundary Wiener diffusion's first passage
#' time, as [wiener()] evaluates it: Navarro and Fuss's (2009) two
#' series combined with a smooth weight, so that the result
#' differentiates exactly on an 'RTMB' tape rather than branching on a
#' parameter. This is the exported seam another extension package builds
#' a joint choice-and-response-time family on, and it is the whole of
#' what this package promises outside itself.
#'
#' @section Tape safety:
#' Every argument may be an 'RTMB' advector except `upper`, which
#' selects a boundary and must be data. The function contains no
#' comparison and no branch on a parameter, so it tapes; the price is
#' that both series are evaluated on every call, which is what makes the
#' density correct over the whole range of normalized times rather than
#' over the part one truncation happens to cover.
#'
#' @section What a decision time of zero or less gives:
#' `-Inf`, not `NaN`. The density is zero there, so the log density is
#' `-Inf`, and that is a value a mixture's log-sum-exp can use. A family
#' calling this is responsible for keeping a fit away from that region;
#' [wiener()] does it with a bounded link on the non-decision time.
#'
#' @param dt Decision time: the response time less the non-decision
#'   time. Not the response time.
#' @param drift Drift rate, on the real line. Positive drift moves the
#'   accumulator toward the upper boundary.
#' @param bs Boundary separation, positive.
#' @param bias Relative start point, in `(0, 1)`. `0.5` starts halfway.
#' @param upper Which boundary the response landed on, `1` for the upper
#'   and `0` for the lower. Data, not a parameter: the reflection that
#'   turns the lower-boundary series into the upper-boundary one is
#'   arithmetic on this rather than a branch.
#'
#' @return A numeric or advector vector, the recycled length of the
#'   arguments.
#'
#' @references
#' Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
#' for first-passage times in Wiener diffusion models. *Journal of
#' Mathematical Psychology* 53, 222-230.
#'
#' @seealso [wiener()], the family this density serves inside this
#'   package.
#'
#' @examples
#' # the density over the upper boundary integrates to the probability
#' # of ever reaching it, which the diffusion has in closed form
#' v <- 1.2; a <- 1.5; w <- 0.3
#' t <- seq(1e-5, 12, length.out = 200000)
#' p_up <- sum(exp(wiener_lpdf(t, v, a, w, 1))) * diff(t)[1]
#' p_lo <- (exp(-2 * v * a) - exp(-2 * v * a * w)) / (exp(-2 * v * a) - 1)
#' c(quadrature = p_up, analytic = 1 - p_lo)
#'
#' # and it differentiates on an RTMB tape, which is the point
#' tp <- RTMB::MakeTape(function(v) sum(wiener_lpdf(c(0.4, 0.9), v,
#'                                                  1.5, 0.3, 1)), 1)
#' tp$jacobian(1.2)
#' @export
wiener_lpdf <- function(dt, drift, bs, bias, upper) {
  if (!inherits(upper, "advector")) {
    u <- as.numeric(upper)
    if (anyNA(u) || !all(u %in% c(0, 1))) {
      stop("wiener_lpdf(upper =) says which boundary each response ",
           "landed on and must be 1 or 0. It is data, so it can be ",
           "checked; the other four arguments are parameters and are ",
           "not", call. = FALSE)
    }
  }
  ddm_lpdf_both(dt, drift, bs, bias, upper)
}
