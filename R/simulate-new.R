# De novo simulation: responses from a bf() specification and supplied
# parameters, no fitted model (the glmmTMB::simulate_new analog; power
# analysis).

#' Simulate responses from a formula and parameters
#'
#' Builds the design from `formula` and `data` exactly as [frm()]
#' would, sets the parameters from `newparams`, and simulates
#' responses. Random effects are redrawn from their covariance for
#' every simulation unless `newparams$b` supplies them.
#'
#' `data` must contain a response column with values that are valid for
#' the family (any dummy values do; they only anchor the design).
#' Inspect the required parameter layout with
#' `frm(formula, data, dry_run = "frame")$par_template`.
#'
#' @param formula A `bf()` formula (with a family attached) or a plain
#'   formula plus `family`.
#' @param data Model data, including a dummy response column.
#' @param family Family, when `formula` does not carry one.
#' @param newparams Named list of parameter vectors matching the
#'   template: `beta` (and `betad`, `theta`, `thetar` as the model
#'   requires); optionally `b` to fix the random effects across
#'   simulations.
#' @param nsim,seed As in [simulate()].
#' @return A data frame with `nsim` columns of simulated responses.
#' @examples
#' # power analysis: simulate from a design with chosen parameters
#' dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)), y = 0)
#' frm(bf(y ~ x + (1 | g)) + gaussian(), dd,
#'     dry_run = "frame")$par_template   # the required layout
#' sims <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
#'                      newparams = list(beta = c(1, 0.5),
#'                                       betad = log(0.7),
#'                                       theta = log(0.5)),
#'                      nsim = 3, seed = 1)
#' head(sims)
#' @export
frm_simulate <- function(formula, data, family = NULL, newparams,
                         nsim = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  bform <- as_bform(formula, family)
  spec <- parse_spec(bform)
  frame <- assemble_frame(spec, data)
  if (length(spec$responses) > 1L) {
    stop("frm_simulate() supports univariate models", call. = FALSE)
  }
  rspec <- spec$responses[[1L]]
  if (is.null(rspec$family$sim)) {
    stop("Family '", rspec$family$family, "' has no simulator yet",
         call. = FALSE)
  }

  est <- frame$par_template
  unknown <- setdiff(names(newparams), c(names(est), "b"))
  if (length(unknown)) {
    stop("Unknown newparams component(s): ",
         paste(unknown, collapse = ", "), " (template has: ",
         paste(names(est), collapse = ", "), ")", call. = FALSE)
  }
  for (nm in names(newparams)) {
    if (length(newparams[[nm]]) != length(est[[nm]])) {
      stop("newparams$", nm, " must have length ", length(est[[nm]]),
           call. = FALSE)
    }
    est[[nm]][] <- newparams[[nm]]
  }
  fixed_b <- "b" %in% names(newparams)

  # a minimal fit-shaped object: eval_dpars and draw_b only touch
  # spec / frame / estimates
  shim <- list(spec = spec, frame = frame, estimates = est)
  av <- frame$aterm_values[[rspec$resp_name]]
  n <- frame$n_obs
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    b_use <- if (is.null(est[["b"]])) {
      NULL
    } else if (fixed_b) {
      est[["b"]]
    } else {
      draw_b(shim)
    }
    dp <- eval_dpars(shim, b = b_use)[[rspec$resp_name]]
    out[[s]] <- rspec$family$sim(dp, av, n)
  }
  names(out) <- paste0("sim_", seq_len(nsim))
  as.data.frame(out)
}
