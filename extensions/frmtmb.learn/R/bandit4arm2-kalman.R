#' A Kalman filter over the restless four-armed bandit
#'
#' The task of Daw and others (2006): four arms whose payoffs drift, so
#' that a subject has to keep exploring. A delta rule with a fixed
#' learning rate is the wrong model for it, because the right learning
#' rate depends on how uncertain the subject currently is about the arm
#' it just pulled. The Kalman filter is that model: it carries a mean
#' AND a variance for each arm, and its gain is the learning rate the
#' uncertainty implies.
#'
#' ```
#' # choice, on the current posterior means
#' P(arm k) <- softmax(tau * mu)
#' # observation, chosen arm only
#' gain     <- v[k] / (v[k] + sigma_o^2)
#' mu[k]    <- mu[k] + gain * (payoff - mu[k])
#' v[k]     <- v[k] * (1 - gain)
#' # diffusion, every arm, ready for the next trial
#' mu       <- lambda * mu + (1 - lambda) * center
#' v        <- lambda^2 * v + sigmaD^2
#' ```
#'
#' hBayesDM calls this `bandit4arm2_kalman_filter` and spells the six
#' parameters `lambda`, `theta`, `beta`, `mu0`, `sigma0` and `sigmaD`.
#'
#' @section What the variance buys, and what it costs:
#' The gain falls as a subject learns, so early trials move the mean
#' further than late ones without any parameter saying so. Nothing else
#' in this package has a learning rate that changes within a subject
#' without a covariate to change it.
#'
#' The cost is identifiability. `sigma0`, `sigmaD` and the fixed
#' observation noise `sigma_o` enter the gain only through ratios, and
#' `tau` trades off against the scale of the payoffs. Fitting all six
#' freely on a short session is not advisable; fix `sigma_o` (it is an
#' argument, not a parameter, for that reason), and consider fixing
#' `mu0` and `center` at the task's own center when the design is known.
#'
#' @section Exploration bonus:
#' `bonus = TRUE` adds Daw and others' exploration bonus: an arm's
#' utility becomes
#'
#' ```
#' util[j] <- tau * (mu[j] + phi * sqrt(s[j]))
#' ```
#'
#' so an arm the subject is UNCERTAIN about is worth more than its
#' posterior mean alone, by an amount `phi` the fit estimates. `phi` is
#' an ordinary distributional parameter, so it takes a formula like any
#' other; `phi = 0` is the plain softmax and is what `bonus = FALSE`,
#' the default, fits.
#'
#' IT IS OFF BY DEFAULT so that every fit made before it existed is
#' unchanged, and the suite pins that by holding `phi` at zero with
#' `bf(phi = 0)` and requiring the two log-likelihoods to agree.
#'
#' WHAT IT TRADES OFF AGAINST IS NOT `tau`, and this help said otherwise
#' before the study ran. The reasoning was that both control how far
#' choice departs from the current best arm. Measured at 30 subjects by
#' 100 trials with `center`, `mu0` and `sigma0` held at the task's own
#' values, 60 replicates: `phi` recovers with a bias of +0.06 against a
#' truth of 1.5, a spread of 0.14 and coverage 0.93, and its estimates
#' correlate with `tau`'s at 0.12, which is nothing. What they correlate
#' with is `sigmaD`, at -0.81, and after the fact that is the obvious
#' pair: `phi` multiplies the posterior standard deviation and `sigmaD`
#' sets how fast that standard deviation grows, so the two scale the
#' same term. Report them together. `tau`'s own interval undercovers a
#' little here, 0.85 against the nominal 0.95.
#'
#' Turn the bonus on when uncertainty-driven exploration is the question
#' the study asks. Leaving `center`, `mu0` and `sigma0` free as well is
#' not advisable on a session this size, for the reason the section
#' above gives.
#'
#' @inheritSection bandit2arm_delta The Laplace caveat
#' @inheritParams bandit2arm_delta
#' @param bonus Add the exploration bonus to the choice rule, giving one
#'   more parameter `phi`. `FALSE`, the default, is the plain softmax
#'   over posterior means and is the model the first release fitted.
#' @param sigma_o The observation noise standard deviation, held FIXED
#'   rather than estimated. Daw and others use 4. It sets the scale the
#'   gain is measured against, and estimating it alongside `sigma0` and
#'   `sigmaD` asks the data to separate three quantities that enter
#'   through two ratios.
#'
#' @return A `frmtmb_family` object, for `frm(family = )`.
#'
#' @references
#' Daw, N. D., O'Doherty, J. P., Dayan, P., Seymour, B. and Dolan, R. J.
#' (2006). Cortical substrates for exploratory decisions in humans.
#' *Nature* 441, 876-879.
#'
#' @seealso [bandit2arm_delta()], [frm_task_design()]
#'
#' @examples
#' d <- frm_task_design("bandit4arm_restless", n_subject = 5, n_trial = 40,
#'                      seed = 6)
#' d$choice <- frm_task_simulate(
#'   bandit4arm2_kalman_filter(subject = id, trial = trial), d,
#'   pars = list(lambda = 0.98, center = 50, tau = 0.15, mu0 = 50,
#'               sigma0 = 10, sigmaD = 3), seed = 6)[[1]]$choice
#' table(d$choice)
#' @export
bandit4arm2_kalman_filter <- function(subject, trial = NULL, sigma_o = 4,
                                      bonus = FALSE) {
  if (!is.numeric(sigma_o) || length(sigma_o) != 1L || is.na(sigma_o) ||
        sigma_o <= 0) {
    stop("bandit4arm2_kalman_filter(sigma_o =) is the observation noise ",
         "standard deviation and must be one positive number; it is held ",
         "fixed rather than estimated", call. = FALSE)
  }
  if (!is.logical(bonus) || length(bonus) != 1L || is.na(bonus)) {
    stop("bandit4arm2_kalman_filter(bonus =) is TRUE or FALSE: whether ",
         "the choice rule adds phi times an arm's posterior standard ",
         "deviation to its value", call. = FALSE)
  }
  vo <- sigma_o^2
  spec <- ln_spec(
    n_option = 4L,
    init = function(ns, d1) {
      one <- rep(1, ns)
      m <- d1[["mu0"]] * one
      v <- d1[["sigma0"]]^2 * one
      list(mu1 = m, mu2 = m, mu3 = m, mu4 = m,
           s1 = v, s2 = v, s3 = v, s4 = v)
    },
    # WITH THE BONUS, the utility is tau * (mu + phi * sigma), which is
    # Daw and others' (2006) "softmax with an exploration bonus": an arm
    # the subject is UNCERTAIN about is worth more than its posterior
    # mean alone, by an amount the fit estimates. sqrt() of a posterior
    # variance is safe on the tape because the variance is a sum of
    # positive terms and cannot reach zero: the diffusion adds sigmaD^2
    # at every trial.
    choice = if (bonus) {
      function(state, d, j) {
        lapply(1:4, function(k) {
          d[["tau"]] * (state[[paste0("mu", k)]] +
                          d[["phi"]] * sqrt(state[[paste0("s", k)]]))
        })
      }
    } else {
      function(state, d, j) {
        lapply(1:4, function(k) d[["tau"]] * state[[paste0("mu", k)]])
      }
    },
    update = function(state, d, ch) {
      lam <- d[["lambda"]]
      vd <- d[["sigmaD"]]^2
      # the payoff the chosen arm returned
      r <- 0
      for (k in 1:4) r <- r + ch[[1L]][[k]] * d[[paste0("payoff", k)]]
      out <- list()
      pe <- 0
      for (k in 1:4) {
        mk <- state[[paste0("mu", k)]]
        vk <- state[[paste0("s", k)]]
        ck <- ch[[1L]][[k]]
        g <- vk / (vk + vo)
        # ck selects the arm that was pulled: an arm not pulled keeps
        # its mean and its variance through the observation step and
        # only diffuses
        out[[paste0("mu", k)]] <- lam * (mk + ck * g * (r - mk)) +
          (1 - lam) * d[["center"]]
        out[[paste0("s", k)]] <- lam^2 * (vk * (1 - ck * g)) + vd
        pe <- pe + ck * (r - mk)
      }
      list(state = out, pe = pe)
    })
  ln_family("bandit4arm2_kalman_filter", substitute(subject),
            substitute(trial),
            # tau is primary, so the main right-hand side of the formula
            # reaches the choice sensitivity. Every other family here
            # puts a learning rate there; this one has no free learning
            # rate, because the gain is derived from the variance, and
            # exploration is what the task was built to measure.
            dpars = c("tau", "lambda", "center", "mu0", "sigma0",
                      "sigmaD", if (bonus) "phi"),
            links = c(list(tau = "log", lambda = "logit",
                           center = "identity", mu0 = "identity",
                           sigma0 = "log", sigmaD = "log"),
                      if (bonus) list(phi = "identity")),
            primary = "tau",
            inits = c(list(tau = function(y, aterms) 0.1,
                           lambda = function(y, aterms) 0.95,
                           center = function(y, aterms) 0,
                           mu0 = function(y, aterms) 0,
                           sigma0 = function(y, aterms) 10,
                           sigmaD = function(y, aterms) 3),
                      if (bonus) list(phi = function(y, aterms) 0)),
            aterms = paste0("payoff", 1:4), spec = spec,
            data_map = stats::setNames(paste0("pay", 1:4),
                                       paste0("payoff", 1:4)),
            constants = list(sigma_o = sigma_o, bonus = bonus))
}
