#' Outcome-representation learning for the Iowa gambling task
#'
#' The model Haines, Vassileva and Ahn (2018) built to fix what
#' [igt_pvl_delta()] cannot see. A prospect-theory learner carries one
#' number per deck, the expected value, so it cannot distinguish a deck
#' that pays a little often from one that pays a lot rarely once the two
#' average out. ORL carries three:
#'
#' * `EV`, the expected value, updated by a delta rule with SEPARATE
#'   rates for a gain and a loss;
#' * `EF`, the expected frequency of a gain, updated by a delta rule on
#'   the SIGN of the outcome, and updated counterfactually on the decks
#'   that were not played;
#' * `PS`, perseverance, which is 1 for the deck just played and decays
#'   toward zero on the others.
#'
#' The three combine linearly into the choice utility,
#'
#' ```
#' util[j] <- EV[j] + betaF * EF[j] + betaP * PS[j]
#' ```
#'
#' with no inverse temperature: the softmax sensitivity is fixed at 1,
#' and the scale of the utilities is carried by the payoffs and by
#' `betaF` and `betaP` instead. That is hBayesDM's parameterization and
#' it is kept, so the payoff SCALE matters here in a way it does not for
#' a family with a free `tau`; see the section below.
#'
#' This is the model hBayesDM calls `igt_orl`.
#'
#' @section How well it is identified, measured:
#' The first release left this family out on the grounds that it is
#' weakly identified. Measured, that claim is half right and the half
#' that is wrong is the half that would have mattered.
#'
#' Every parameter RECOVERS. At 30 subjects by 100 trials with a random
#' intercept on `Arew`, 60 replicates, every bias is at or inside its
#' Monte Carlo error and coverage runs 0.92 to 0.98; 60 of 60 fits gave
#' a usable interval. So a point estimate from this family is not the
#' hazard a variance component from a twenty-trial session is.
#'
#' What is weak is the SEPARATION of the three choice-rule parameters,
#' and it shows up in the correlations of the estimates rather than in
#' their biases. Across replicates `betaF` and `betaP` correlate at
#' -0.77, `k` with `betaP` at -0.53 and `k` with `betaF` at +0.45. The
#' two learning rates are comparatively clean (`Arew` with `Apun`,
#' 0.38) and neither correlates with the choice-rule three above 0.21.
#'
#' The practical consequence is specific. A study comparing groups on
#' `Arew` or `Apun` is on solid ground. A study comparing them on
#' `betaF` alone is at risk of attributing to outcome frequency what
#' belongs to perseverance, because the two are estimated against each
#' other; report them together, or test them jointly. `?frmtmb.learn`
#' carries the whole table.
#'
#' @section The payoff scale:
#' hBayesDM divides the Iowa gambling task's payoffs by 100 before
#' fitting, and the reason is not cosmetic: with the softmax sensitivity
#' fixed at 1, doubling the payoffs doubles every utility and makes the
#' model twice as deterministic. `frm_task_design("igt")` already
#' returns payoffs on that scale. Data brought in raw dollars must be
#' divided the same way, or `betaF` and `betaP` are estimated against a
#' `EV` a hundred times larger than they were built for.
#'
#' @section What differs from hBayesDM:
#' One transform. hBayesDM puts `K` on `(0, 5)` through a probit and
#' decays perseverance by `3^K`; this family estimates `k` on
#' `(0, Inf)` with a log link and decays by the same `3^k`. So the two
#' `k` values are the same number on the same scale, and only the bound
#' differs: a fit here can report a decay faster than hBayesDM's prior
#' allows, and a fit that wants to must be reported as such rather than
#' silently pinned at 5. Everything else is a rename: `Arew` and `Apun`
#' keep their names, and `betaF` and `betaP` keep theirs.
#'
#' @section The data a trial carries:
#' `payoff(pay1, pay2, pay3, pay4)`, one column per deck, holding what
#' each deck WOULD have paid on that trial. The likelihood reads the
#' played deck's entry; the other three make the simulator coherent, and
#' the counterfactual frequency update reads only the SIGN the played
#' deck produced, not the others' values.
#'
#' @param subject The column separating one learner's trial sequence
#'   from the next, given unquoted.
#' @param trial The column giving trial order within a subject, given
#'   unquoted. `NULL` uses the order the rows appear in.
#'
#' @return A `frmtmb_family` object, for `frm(family = )`.
#'
#' @references
#' Haines, N., Vassileva, J. and Ahn, W.-Y. (2018). The
#' outcome-representation learning model: a novel reinforcement learning
#' model of the Iowa gambling task. *Cognitive Science* 42, 2534-2561.
#'
#' @seealso [igt_pvl_delta()] for the prospect-theory model of the same
#'   task, [frm_value_trace()] for the three value stores per trial.
#'
#' @examples
#' d <- frm_task_design("igt", n_subject = 6, n_trial = 60, seed = 8)
#' d$choice <- frm_task_simulate(
#'   igt_orl(subject = id, trial = trial), d,
#'   pars = list(Arew = 0.3, Apun = 0.1, k = 0.5, betaF = 1,
#'               betaP = 1), seed = 8)[[1]]$choice
#' table(d$choice)
#' @export
igt_orl <- function(subject, trial = NULL) {
  spec <- ln_spec(
    n_option = 4L,
    init = function(ns, d1) {
      z <- rep(0, ns)
      list(ev1 = z, ev2 = z, ev3 = z, ev4 = z,
           ef1 = z, ef2 = z, ef3 = z, ef4 = z,
           ps1 = z, ps2 = z, ps3 = z, ps4 = z)
    },
    # no inverse temperature: the sensitivity is fixed at 1 and the
    # scale lives in the payoffs and in the two weights
    choice = function(state, d, j) {
      lapply(1:4, function(k) {
        state[[paste0("ev", k)]] +
          d[["betaF"]] * state[[paste0("ef", k)]] +
          d[["betaP"]] * state[[paste0("ps", k)]]
      })
    },
    update = function(state, d, ch) {
      # The outcome the played deck returned, and its sign. Both are
      # DATA: the payoffs are addition-term columns and the indicator
      # that picks one is the observed choice, so the branch on the sign
      # is a multiplication by a column of zeros and ones and reaches
      # the tape as arithmetic rather than as a comparison.
      x <- 0
      for (k in 1:4) x <- x + ch[[1L]][[k]] * d[[paste0("payoff", k)]]
      pos <- as.numeric(x >= 0)
      s <- sign(x)
      # the rate that moves the played deck, and the one that moves the
      # other three counterfactually: they SWAP on a loss, which is the
      # asymmetry the model is named for
      a_play <- pos * d[["Arew"]] + (1 - pos) * d[["Apun"]]
      a_fic <- pos * d[["Apun"]] + (1 - pos) * d[["Arew"]]
      # hBayesDM's K reaches the recursion as a decay of 3^K, through
      # PS / (1 + K_tr) with K_tr = 3^K - 1. Written as the division it
      # is, so the reader is not asked to cancel the two ones.
      dec <- exp(d[["k"]] * log(3))
      out <- list()
      pe <- 0
      for (k in 1:4) {
        ck <- ch[[1L]][[k]]
        ev <- state[[paste0("ev", k)]]
        ef <- state[[paste0("ef", k)]]
        ps <- state[[paste0("ps", k)]]
        out[[paste0("ev", k)]] <- ev + ck * a_play * (x - ev)
        # played: toward the sign of what it paid. Unplayed: toward
        # minus a third of that sign, which is what spreads one deck's
        # news over the three that were not seen
        out[[paste0("ef", k)]] <- ef + ck * a_play * (s - ef) +
          (1 - ck) * a_fic * (-s / 3 - ef)
        out[[paste0("ps", k)]] <- ck + (1 - ck) * ps / dec
        pe <- pe + ck * (x - ev)
      }
      list(state = out, pe = pe, outcome = x)
    })
  ln_family("igt_orl", substitute(subject), substitute(trial),
            dpars = c("Arew", "Apun", "k", "betaF", "betaP"),
            links = list(Arew = "logit", Apun = "logit", k = "log",
                         betaF = "identity", betaP = "identity"),
            # Arew is primary, so the main right-hand side reaches the
            # learning rate for gains, which is the parameter the
            # clinical literature compares between groups.
            primary = "Arew",
            inits = list(Arew = function(y, aterms) 0.3,
                         Apun = function(y, aterms) 0.1,
                         k = function(y, aterms) 0.5,
                         betaF = function(y, aterms) 0,
                         betaP = function(y, aterms) 0),
            aterms = paste0("payoff", 1:4), spec = spec,
            data_map = stats::setNames(paste0("pay", 1:4),
                                       paste0("payoff", 1:4)))
}
