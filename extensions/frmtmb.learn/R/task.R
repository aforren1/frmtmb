#' Build a trial-level design for a learning task
#'
#' One row per subject and trial, with the payoff schedule of EVERY
#' option fixed in advance and a placeholder response. Fixing the whole
#' schedule ahead of the choices is what makes a draw from these models
#' coherent: a simulated subject that takes the arm the real one did not
#' still needs to be paid.
#'
#' The column names are canonical, because [frm_task_simulate()] reads
#' them: `id`, `trial`, `choice`, and `pay1` to `pay2` or `pay4`. A fit
#' to your own data needs none of this; the family reads whatever
#' columns your formula names.
#'
#' @param task One of `"bandit2arm"` (a stationary two-armed bandit),
#'   `"reversal"` (the same with the contingency reversed at the
#'   halfway point, and an `after_reversal` factor to model it with),
#'   `"bandit4arm_restless"` (four arms whose payoffs follow the
#'   decaying Gaussian random walk of Daw et al. 2006), `"igt"` (the
#'   Iowa gambling task's four decks) or `"twostep"` (the two-stage
#'   Markov task of Daw et al. 2011).
#' @param n_subject,n_trial Subjects, and trials per subject.
#' @param p Two payoff probabilities, for the two-armed tasks.
#' @param seed Passed to [set.seed()] when not `NULL`.
#' @param ... Task-specific settings: `decay` and `center` and `sd_walk`
#'   and `sd_obs` for the restless bandit, `p_reward` for the two-step
#'   task's stage-two options.
#'
#' @return A data frame, ready for [frm_task_simulate()].
#'
#' @seealso [frm_task_simulate()], [bandit2arm_delta()]
#'
#' @examples
#' d <- frm_task_design("reversal", n_subject = 4, n_trial = 20, seed = 1)
#' head(d)
#' table(d$after_reversal)
#' @export
frm_task_design <- function(task = c("bandit2arm", "reversal",
                                     "bandit4arm_restless", "igt",
                                     "twostep"),
                            n_subject = 30L, n_trial = 100L,
                            p = c(0.7, 0.3), seed = NULL, ...) {
  task <- match.arg(task)
  if (!is.null(seed)) set.seed(seed)
  dots <- list(...)
  n_subject <- as.integer(n_subject)
  n_trial <- as.integer(n_trial)
  if (n_subject < 1L || n_trial < 1L) {
    stop("frm_task_design(): n_subject and n_trial must both be at ",
         "least 1", call. = FALSE)
  }
  d <- expand.grid(trial = seq_len(n_trial), id = seq_len(n_subject))
  # expand.grid() leaves an out.attrs attribute that str() prints at
  # length and nothing reads
  attr(d, "out.attrs") <- NULL
  d <- d[, c("id", "trial")]
  n <- nrow(d)
  d$id <- factor(d$id)
  d$choice <- 1L
  if (task %in% c("bandit2arm", "reversal")) {
    pr <- matrix(rep(p, each = n), n, 2L)
    if (task == "reversal") {
      after <- d$trial > (n_trial %/% 2L)
      pr[after, ] <- pr[after, 2:1]
      d$after_reversal <- factor(after, c(FALSE, TRUE), c("before", "after"))
    }
    d$pay1 <- stats::rbinom(n, 1L, pr[, 1L])
    d$pay2 <- stats::rbinom(n, 1L, pr[, 2L])
  } else if (task == "bandit4arm_restless") {
    # Daw et al. (2006): each arm's mean decays toward a common center
    # and diffuses, and the payoff is that mean plus observation noise.
    decay <- dots[["decay"]] %||% 0.9836
    center <- dots[["center"]] %||% 50
    sd_walk <- dots[["sd_walk"]] %||% 2.8
    sd_obs <- dots[["sd_obs"]] %||% 4
    for (k in 1:4) d[[paste0("pay", k)]] <- NA_real_
    for (s in levels(d$id)) {
      rows <- which(d$id == s)
      rows <- rows[order(d$trial[rows])]
      mu <- stats::runif(4L, center - 20, center + 20)
      for (t in seq_along(rows)) {
        for (k in 1:4) d[[paste0("pay", k)]][rows[t]] <-
          stats::rnorm(1L, mu[k], sd_obs)
        mu <- decay * mu + (1 - decay) * center + stats::rnorm(4L, 0, sd_walk)
      }
    }
  } else if (task == "igt") {
    # The four decks of the original task, in units of 100 as the
    # PVL-delta literature scales them: two bad decks with large gains
    # and rare larger losses, two good decks with small gains.
    gain <- c(1, 1, 0.5, 0.5)
    loss <- c(-2.5, -12.5, -0.5, -2.5)
    lossp <- c(0.5, 0.1, 0.5, 0.1)
    for (k in 1:4) {
      d[[paste0("pay", k)]] <- gain[k] +
        loss[k] * stats::rbinom(n, 1L, lossp[k])
    }
  } else {
    # two-step: the four stage-two options drift as independent bounded
    # random walks, which is how the task keeps model-free and
    # model-based learners apart
    pw <- dots[["p_reward"]] %||% c(0.6, 0.4, 0.4, 0.6)
    for (k in 1:4) d[[paste0("pay", k)]] <- NA_real_
    for (s in levels(d$id)) {
      rows <- which(d$id == s)
      rows <- rows[order(d$trial[rows])]
      q <- pw
      for (t in seq_along(rows)) {
        for (k in 1:4) d[[paste0("pay", k)]][rows[t]] <-
          stats::rbinom(1L, 1L, q[k])
        q <- pmin(0.75, pmax(0.25, q + stats::rnorm(4L, 0, 0.025)))
      }
    }
    d$state2 <- 1L
    d$choice2 <- 1L
  }
  rownames(d) <- NULL
  d
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Draw datasets from a learning family's generative process
#'
#' Walks each subject forward through the design, drawing a choice from
#' the current value store and learning from the payoff that choice
#' earns. It is the same recursion the likelihood tapes, at
#' `mode = "simulate"`, so a draw and the density that scores it cannot
#' drift apart.
#'
#' This route takes parameters directly, per subject or per row, rather
#' than through a formula and a random-effect covariance. That is what a
#' recovery study wants, and it makes the simulator independent of the
#' fitting machinery, so agreement between a fit and its truth is
#' evidence about both. [frmtmb::frm_simulate()] is the other route and
#' goes through the fitted grammar; the two agree by construction for
#' every family except [ts_par7()], which [frmtmb::frm_simulate()]
#' refuses because one trial's draw is three numbers.
#'
#' @param family A family from this package, such as
#'   [bandit2arm_delta()].
#' @param data A design from [frm_task_design()], or any data frame with
#'   the same canonical columns.
#' @param pars Named list of parameters on their NATURAL scales (a
#'   learning rate in `(0, 1)`, not its logit). Each element is one
#'   value, recycled everywhere; or one value per subject, in the order
#'   `levels(factor(data$id))` gives; or one value per ROW. The per-row
#'   form is what a covariate that varies WITHIN a subject needs, which
#'   is what a reversal is: see the example.
#' @param nsim Datasets to draw.
#' @param response Name of the column the drawn response goes into.
#'   `NULL`, the default, uses the family's own: `choice` for a family
#'   whose response is the option taken, `rt` for [rlddm()], whose
#'   response is the response time and whose drawn choice is written
#'   beside it.
#' @param seed Passed to [set.seed()] when not `NULL`.
#'
#' @return A list of `nsim` data frames, each `data` with the drawn
#'   columns replaced.
#'
#' @seealso [frm_task_design()], [frmtmb::frm_simulate()]
#'
#' @examples
#' d <- frm_task_design("bandit2arm", n_subject = 5, n_trial = 30,
#'                      seed = 1)
#' sims <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
#'                           d, pars = list(alpha = 0.4, tau = 3),
#'                           nsim = 2, seed = 1)
#' table(sims[[1]]$choice)
#'
#' # a learning rate that changes within a subject, one value per row
#' r <- frm_task_design("reversal", n_subject = 5, n_trial = 30, seed = 2)
#' a <- ifelse(r$after_reversal == "after", 0.6, 0.2)
#' rs <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
#'                         r, pars = list(alpha = a, tau = 3), seed = 2)
#' table(rs[[1]]$choice)
#' @export
frm_task_simulate <- function(family, data, pars, nsim = 1L,
                              response = NULL, seed = NULL) {
  lrn <- family[["learn"]]
  if (is.null(lrn)) {
    stop("frm_task_simulate() takes a family from frmtmb.learn, and ",
         "this one carries no learning recursion", call. = FALSE)
  }
  nm <- family[["family"]]
  n <- nrow(data)
  gv <- factor(eval(lrn[["subject_expr"]], data, parent.frame()))
  tv <- if (is.null(lrn[["trial_expr"]])) {
    out <- integer(n)
    for (g in split(seq_len(n), gv)) out[g] <- seq_along(g)
    out
  } else {
    as.numeric(eval(lrn[["trial_expr"]], data, parent.frame()))
  }
  block <- ln_pack(gv, tv, n, nm)
  miss <- setdiff(lrn[["dpars"]], names(pars))
  if (length(miss)) {
    stop("frm_task_simulate(): ", nm, "() has no value for ",
         paste(miss, collapse = ", "), ". Every parameter of the family ",
         "needs one, on its natural scale", call. = FALSE)
  }
  si <- as.integer(gv)
  cd <- list()
  for (dp in lrn[["dpars"]]) {
    v <- as.numeric(pars[[dp]])
    # Three lengths, and the per-ROW one is not a convenience: a
    # covariate that varies WITHIN a subject, which is what a reversal
    # is, gives a parameter that changes from trial to trial. A
    # simulator that took only per-subject values could not draw the
    # data this package's own worked example fits.
    if (length(v) == 1L) {
      cd[[dp]] <- v
    } else if (length(v) == block[["n_subj"]] && block[["n_subj"]] != n) {
      cd[[dp]] <- v[si]
    } else if (length(v) == n) {
      cd[[dp]] <- v
    } else {
      stop("frm_task_simulate(): '", dp, "' has ", length(v), " values ",
           "for ", block[["n_subj"]], " subjects and ", n, " rows. Give ",
           "one value, one per subject in the order levels(factor(id)) ",
           "gives, or one per row", call. = FALSE)
    }
  }
  dm <- lrn[["data_map"]]
  need <- setdiff(unname(dm), names(data))
  if (length(need)) {
    stop("frm_task_simulate(): the design has no column '", need[1L],
         "'. frm_task_design() names the payoff columns pay1, pay2 and ",
         "so on, and this route reads them by those names", call. = FALSE)
  }
  for (k in names(dm)) cd[[k]] <- as.numeric(data[[dm[[k]]]])
  # A family whose choice rule is a DENSITY has a response that is not
  # an option code: rlddm()'s is a response time, which is continuous
  # and belongs in a column called rt rather than choice. Both facts
  # follow from the one property, so both are read off it rather than
  # declared twice.
  # A family whose DRAW needs a package the likelihood does not, checked
  # once here rather than once per trial inside the walk. rlddm() is the
  # one: its joint draw of a boundary and a time is exact only through
  # RWiener, and refusing up front is better than the alternative,
  # which is a draw that succeeds on most rows and errors on the rest.
  need <- lrn[["sim_needs"]]
  if (!is.null(need) && !requireNamespace(need, quietly = TRUE)) {
    stop("frm_task_simulate(): drawing from ", nm, "() needs the '",
         need, "' package, which is not installed. It draws a boundary ",
         "and a response time jointly from the diffusion, and the ",
         "route that does so exactly is the one '", need, "' supplies. ",
         "Fitting needs nothing extra: the likelihood is exact without ",
         "it", call. = FALSE)
  }
  code_resp <- is.null(lrn[["spec"]][["logp"]])
  if (is.null(response)) {
    response <- if (code_resp) "choice" else lrn[["sim_response"]]
  }
  if (!is.null(seed)) set.seed(seed)
  lapply(seq_len(nsim), function(s) {
    out <- ln_recurse(block, cd, rep(1, n), lrn[["spec"]], "simulate")
    d <- data
    d[[response]] <- if (code_resp) as.integer(out[["y"]]) else out[["y"]]
    for (nmc in names(out[["cols"]])) {
      col <- lrn[["sim_out"]][[nmc]]
      if (!is.null(col)) d[[col]] <- as.integer(out[["cols"]][[nmc]])
    }
    d
  })
}
