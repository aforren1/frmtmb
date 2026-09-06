#' Prospect-valence learning with a delta rule, for the Iowa gambling
#' task
#'
#' Four decks, each paying a gain on every trial and an occasional
#' larger loss. What makes this model different from the bandits is that
#' the OUTCOME is transformed before it is learned from: money is
#' converted to subjective value by a prospect-theory utility with a
#' compressive exponent and a loss-aversion multiplier, and the delta
#' rule then runs on that.
#'
#' ```
#' u <- if (x >= 0)  abs(x)^shape
#'      else        -lambda * abs(x)^shape
#' Q[chosen] <- Q[chosen] + alpha * (u - Q[chosen])
#' P(deck k) <- softmax(tau * Q)
#' ```
#'
#' @section Names, and one trap:
#' hBayesDM calls this `igt_pvl_delta` and spells the parameters `A`,
#' `alpha`, `lambda` and `cons`. The map is
#'
#' | this family | hBayesDM | what it is |
#' |---|---|---|
#' | `alpha` | `A` | the learning rate, on (0, 1) |
#' | `shape` | `alpha` | the utility exponent, on (0, 1); see below |
#' | `lambda` | `lambda` | loss aversion, positive |
#' | `tau` | `3^cons - 1` | softmax sensitivity, positive |
#'
#' **`alpha` means different things in the two packages.** It is the
#' learning rate here, as it is in every other family in this package,
#' and it is the utility exponent in hBayesDM. Reading an hBayesDM
#' `alpha` into this family's `alpha` fits a different model that
#' converges quietly.
#'
#' **`shape` reaches a NARROWER range than hBayesDM's `alpha`.** The
#' logit link confines it to (0, 1); hBayesDM bounds the same quantity
#' on (0, 2), through `Phi_approx(...) * 2`. An hBayesDM fit whose
#' utility exponent came out above 1, which is a convex rather than a
#' compressive utility, is not representable here. Every other bound
#' difference between the two runs the other way, this package being the
#' more permissive, so this is the one to check before porting.
#'
#' `tau` is the fourth difference and it is a reparameterization rather
#' than a rename: hBayesDM estimates a consistency parameter `cons` on
#' `(0, 5)` and uses `3^cons - 1` as the sensitivity, while this family
#' estimates the sensitivity itself on a log link. The two are the same
#' model at `tau = 3^cons - 1`, and `tau` is then the same quantity, on
#' the same scale, that it is in every other family here.
#'
#' @section Why PVL-delta and not ORL:
#' ORL (Haines, Vassileva and Ahn 2018) is the other Iowa gambling model
#' worth having, and it is not here. PVL-delta was chosen because its
#' four parameters recover at a realistic scale, which is the standard
#' this package holds a family to, and because its outcome transform
#' exercises a part of the engine nothing else does: a nonlinear
#' function of DATA raised to an estimated exponent. ORL carries two
#' further value stores and two frequency-weighting parameters that are
#' known to be weakly identified, so shipping it would mean shipping a
#' recovery table that says so. The engine takes it in about forty
#' lines when someone wants it.
#'
#' @section Scale:
#' The exponent makes this family sensitive to the units of the payoffs.
#' `frm_task_design("igt")` returns them in units of 100, which is the
#' scaling the PVL literature uses; passing raw currency changes what
#' `shape` and `tau` mean and will usually make the fit worse
#' conditioned.
#'
#' @inheritSection bandit2arm_delta The Laplace caveat
#' @inheritParams bandit2arm_delta
#'
#' @return A `frmtmb_family` object, for `frm(family = )`.
#'
#' @references
#' Ahn, W.-Y., Busemeyer, J. R., Wagenmakers, E.-J. and Stout, J. C.
#' (2008). Comparison of decision learning models using the generalization
#' criterion method. *Cognitive Science* 32, 1376-1402.
#'
#' @seealso [bandit4arm2_kalman_filter()], [frm_task_design()]
#'
#' @examples
#' d <- frm_task_design("igt", n_subject = 6, n_trial = 50, seed = 7)
#' d$choice <- frm_task_simulate(
#'   igt_pvl_delta(subject = id, trial = trial), d,
#'   pars = list(alpha = 0.3, shape = 0.4, lambda = 1.5, tau = 1),
#'   seed = 7)[[1]]$choice
#' # the good decks are 3 and 4
#' table(d$choice)
#' @export
igt_pvl_delta <- function(subject, trial = NULL) {
  spec <- ln_spec(
    n_option = 4L,
    init = function(ns, d1) {
      z <- rep(0, ns)
      list(q1 = z, q2 = z, q3 = z, q4 = z)
    },
    choice = function(state, d, j) {
      lapply(1:4, function(k) d[["tau"]] * state[[paste0("q", k)]])
    },
    update = function(state, d, ch) {
      # the outcome the chosen deck returned, which is DATA, so its
      # sign and its magnitude are both data and only the exponent and
      # the multiplier are estimated
      x <- 0
      for (k in 1:4) x <- x + ch[[1L]][[k]] * d[[paste0("payoff", k)]]
      w <- abs(x)
      nz <- as.numeric(w > 0)
      # log(0) is -Inf and 0 * -Inf is NaN, so a zero outcome is raised
      # to the power of 1 and then multiplied out again
      mag <- exp(d[["shape"]] * log(w + (1 - nz))) * nz
      u <- (as.numeric(x >= 0) -
              (1 - as.numeric(x >= 0)) * d[["lambda"]]) * mag
      out <- list()
      pe <- 0
      for (k in 1:4) {
        qk <- state[[paste0("q", k)]]
        ck <- ch[[1L]][[k]]
        out[[paste0("q", k)]] <- qk + ck * d[["alpha"]] * (u - qk)
        pe <- pe + ck * (u - qk)
      }
      list(state = out, pe = pe, utility = u)
    })
  ln_family("igt_pvl_delta", substitute(subject), substitute(trial),
            dpars = c("alpha", "shape", "lambda", "tau"),
            links = list(alpha = "logit", shape = "logit",
                         lambda = "log", tau = "log"),
            primary = "alpha",
            inits = list(alpha = function(y, aterms) 0.3,
                         shape = function(y, aterms) 0.5,
                         lambda = function(y, aterms) 1,
                         tau = function(y, aterms) 1),
            aterms = paste0("payoff", 1:4), spec = spec,
            data_map = stats::setNames(paste0("pay", 1:4),
                                       paste0("payoff", 1:4)))
}
