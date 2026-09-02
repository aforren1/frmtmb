# Ordinary differential equation dynamics inside a nonlinear predictor.
#
# The whole feature is one exported helper. bf(nl = TRUE) already
# evaluates an arbitrary R body on the AD tape, and RTMBode::ode() is
# already differentiable through both its initial states and its
# parameters, so nothing in the parser, the frame, or the objective has
# to learn about ODEs. What a user cannot reasonably be asked to write
# by hand is the per-group solve loop and its failure modes, which is
# what frm_ode() owns.
#
# The one hard constraint (dev/ode-feasibility.md, section 7): solve one
# small system per group, never one stacked system over all groups. The
# Laplace inner problem's second-order path through the adjoint ODE node
# returns NaN gradients above roughly 8 states with lsoda, and crashes
# the process at 32 states with some integrators. Per-group compartment
# models live well below that ceiling; a stacked system does not. No
# future optimization may collapse the loop.

# Adaptive integrators only. A fixed-step integrator run at deSolve's
# default step count does not solve the system to the tolerance the
# likelihood is defined at, so it returns a *different* likelihood, not
# a slightly noisier one: rk4 reported -63.93 where every adaptive
# method reported -60.46 on the same probe data.
ode_adaptive_methods <- c(
  "lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk", "radau",
  "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "ode45", "ode23"
)

# deSolve says "integration was not successful" through a warning and
# then returns a full-length matrix of finite garbage, so a solve that
# gave up is otherwise indistinguishable from one that worked.
ode_giveup_pattern <- paste(
  "not successful", "excessive amount of work", "maxsteps",
  "Returning early", "singular", "step size", "convergence",
  sep = "|"
)

# Above this many states in a single system, the Laplace gradient
# through the adjoint node starts returning NaN (integrator-dependent).
ode_state_warn <- 8L

# One-shot warning bookkeeping: the nl body is taped, so the R code
# below runs once or twice per fit, but re-fitting in the same session
# should not repeat a structural warning the user has already read.
ode_warned <- new.env(parent = emptyenv())

# The groups whose solve failed on the most recent frm_ode() call, so
# that a penalty written into a prediction can be read back after the
# fact. Reset by every call; see frm_ode_failures().
ode_failure_log <- new.env(parent = emptyenv())

ode_warn_once <- function(key, ...) {
  if (!is.null(ode_warned[[key]])) return(invisible(NULL))
  ode_warned[[key]] <- TRUE
  warning(..., call. = FALSE)
}

#' Require the optional solver backend.
#'
#' @noRd
ode_has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

ode_require_backend <- function() {
  pkgs <- c("RTMBode", "deSolve")
  miss <- pkgs[!vapply(pkgs, ode_has_pkg, TRUE)]
  if (!length(miss)) return(invisible(TRUE))
  stop("frm_ode() needs the ", paste(miss, collapse = " and "),
       " package", if (length(miss) > 1L) "s" else "",
       ", which ", if (length(miss) > 1L) "are" else "is",
       " not installed. RTMBode is not on CRAN; install it with\n",
       "  install.packages(\"RTMBode\", repos = c(\n",
       "    \"https://kaskr.r-universe.dev\",\n",
       "    \"https://cloud.r-project.org\"))",
       call. = FALSE)
}

#' Normalize `init` / `parms` to a list of columns.
#'
#' Columns, not rows: each element is either one value per observation
#' (a linear predictor or a data column, read off the group's first row)
#' or a single value shared by every group. A bare vector is one column,
#' never a row of per-state values - the two readings are impossible to
#' tell apart when the number of states happens to equal the number of
#' observations, so the ambiguous spelling is refused instead of
#' guessed.
#'
#' @noRd
ode_columns <- function(x, n_obs, arg) {
  if (is.list(x) && !inherits(x, "advector")) {
    cols <- as.list(x)
  } else if (!is.null(dim(x))) {
    d <- dim(x)
    if (d[1L] != n_obs) {
      stop("`", arg, "` is a matrix with ", d[1L], " rows but there are ",
           n_obs, " observations", call. = FALSE)
    }
    cols <- lapply(seq_len(d[2L]), function(j) x[, j])
  } else {
    cols <- list(x)
  }
  if (!length(cols)) {
    stop("`", arg, "` is empty; it needs at least one column",
         call. = FALSE)
  }
  for (j in seq_along(cols)) {
    len <- length(cols[[j]])
    if (len != 1L && len != n_obs) {
      stop("`", arg, "` column ", j, " has length ", len,
           "; it must be length ", n_obs,
           " (one value per observation) or length 1 (shared by every ",
           "group). To give several values, pass a list: ", arg,
           " = list(a, b, ...)", call. = FALSE)
    }
  }
  cols
}

#' Within-group constancy of the columns that arrive as plain data.
#'
#' A dynamics parameter that varies inside a solve group is a different
#' model from the one frm_ode() can solve, and the fit that results is
#' silently wrong rather than loud (the coefficient stays at its start
#' value and the Hessian goes indefinite). Columns carrying AD values
#' cannot be inspected here - RTMB refuses comparison on AD types, by
#' design - so those are checked structurally at frame assembly by
#' `check_ode_constancy()`. Plain numeric columns are checked here,
#' where the actual values are available.
#'
#' @noRd
ode_check_constant <- function(cols, groups, arg, labels) {
  for (j in seq_along(cols)) {
    v <- cols[[j]]
    if (inherits(v, "advector") || length(v) == 1L || !is.numeric(v)) {
      next
    }
    for (g in seq_along(groups)) {
      idx <- groups[[g]]
      if (length(idx) < 2L) next
      rng <- range(v[idx])
      if (rng[2L] - rng[1L] >
            1e-8 * max(1, max(abs(rng)))) {
        stop("`", arg, "` column ", j, " is not constant within group '",
             labels[[g]], "' (values ", format(rng[1L]), " to ",
             format(rng[2L]), "). frm_ode() solves one system per group ",
             "and reads each dynamics input off the group's first row, ",
             "so a within-group covariate cannot enter the likelihood.",
             call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

#' Solve an ODE once per group inside a nonlinear predictor
#'
#' Evaluates a system of ordinary differential equations separately for
#' each group of rows and returns the solution aligned with the rows of
#' the data, so that a compartment model can be written directly in the
#' body of a `bf(..., nl = TRUE)` formula. The dynamics parameters and
#' the initial states are ordinary nonlinear parameters, which means
#' they take fixed effects, random effects and covariates like any other
#' linear predictor, and the Laplace approximation is exact through the
#' solver's adjoint.
#'
#' @details
#' # What the group is
#'
#' `group` names the unit that owns one ODE system: a subject in a
#' population pharmacokinetic model, a reactor, a patient. Every row of
#' one group is one observation of that group's trajectory, at the time
#' given by `times`. One solve is performed per group, over that group's
#' own sorted times, and the results are scattered back into the input
#' row order. Ragged designs, unsorted rows, repeated times and an
#' observation at `t0` itself are all fine.
#'
#' Groups are never stacked into one large system. That is a hard
#' constraint, not an implementation detail: the second-order derivative
#' path the Laplace approximation needs returns `NaN` above roughly
#' eight states in a single system, so a stacked solve gives silently
#' wrong or missing gradients. `frm_ode()` warns when one system alone
#' exceeds that many states.
#'
#' # What is constant within a group
#'
#' `init` and `parms` are read off each group's **first row**. A
#' dynamics parameter must therefore not vary inside a group: a
#' covariate that changes between an early and a late observation of the
#' same subject describes a model this helper cannot solve. Such a
#' covariate is refused, by name, rather than silently ignored.
#' Covariates that are constant within a group (a subject's weight, dose
#' or treatment arm) are the intended case and are unrestricted.
#'
#' # Failed solves
#'
#' `on_error` and `penalty` reach less than they look like they do, and
#' the difference matters:
#'
#' - **While fitting**, the body is evaluated on the
#'   automatic-differentiation tape. `RTMBode::ode()` returns an
#'   `advector` and does not raise an R error when the trajectory goes
#'   bad, and a non-finite AD value cannot be tested for (RTMB refuses
#'   comparison on AD types). So a diverging region of the parameter
#'   space is **not** caught here. It surfaces instead as the
#'   optimizer's own `NA/NaN function evaluation` or `NA/NaN gradient
#'   evaluation` warning, and `on_error = "error"` will not name the
#'   group. The two failures that are still caught on the tape are an
#'   integrator that gives up and returns fewer time points than were
#'   asked for, and any R error raised by `dynamics` itself.
#' - **Everywhere else** the evaluation is ordinary numeric arithmetic,
#'   and every check applies: [predict()], [simulate()],
#'   [residuals()], a direct call, and a fit whose `init` and `parms`
#'   contain no estimated parameter (which is evaluated once, numerically,
#'   as the tape is built). Here `on_error = "penalize"` writes `penalty`
#'   into the failed group's rows, and `on_error = "error"` stops and
#'   names the group.
#'
#' A penalty is never written silently: `frm_ode()` warns, naming the
#' groups, and [frm_ode_failures()] reports them afterwards. Treat those
#' rows as missing, not as predictions.
#'
#' If a fit reports `NA/NaN gradient evaluation`, run `diagnose()` on it,
#' then call `frm_ode()` directly at the suspect parameter values with
#' `on_error = "error"`: numerically it will name the group that cannot
#' be solved.
#'
#' Solver warnings from \pkg{deSolve} ("corrector convergence failed
#' repeatedly", "exceeded maxsteps") during a fit come from the
#' optimizer's probing steps and are usually not fatal. Judge the fit by
#' the gradient at the optimum, not by whether the solver complained.
#'
#' # Boundaries
#'
#' Dosing event tables (`evid`/`amt` records, repeated doses,
#' infusions) are out of scope. Only models driven by their initial
#' conditions are supported. For event-driven population
#' pharmacokinetics use a dedicated tool such as `nlmixr2`.
#'
#' `predict(se.fit = TRUE)` is not available for a nonlinear predictor,
#' including one containing `frm_ode()`; request a nonlinear parameter
#' with `predict(dpar = )` instead.
#'
#' @param dynamics A function `function(t, y, parms)` giving the
#'   derivatives, following the \pkg{deSolve} convention: `t` is the
#'   scalar time, `y` the state vector, `parms` the parameter vector,
#'   and the return value is `list(dydt)`. A bare derivative vector is
#'   also accepted. Index `y` and `parms` by position; use
#'   `"c" <- RTMB::ADoverload("c")` inside the function so that `c()`
#'   keeps the automatic-differentiation class.
#' @param init Initial states, one column per state: a list, a matrix,
#'   or a single vector for a one-state system. Each column is either
#'   one value per observation (constant within group) or one value
#'   shared by every group. The number of columns sets the number of
#'   states.
#' @param times Observation times, one per row of the data.
#' @param parms Dynamics parameters, one column per parameter, in the
#'   order `dynamics` expects them. Same shape rules as `init`.
#' @param group Grouping vector naming the unit that owns one system,
#'   one value per row. `NULL` treats the whole data as one group.
#' @param output Which states to return: `NULL` (the default) returns
#'   every state, an integer or character vector selects some. A single
#'   selected state is returned as a vector, otherwise a matrix with one
#'   column per selected state. Character selection requires `states`.
#' @param states Optional state names, one per column of `init`.
#' @param t0 Initial time, a scalar or one value per row (constant
#'   within group). Every observation time must be at or after it.
#' @param method Integrator, passed to [deSolve::ode()]. Must be
#'   adaptive; fixed-step integrators such as `"rk4"` and `"euler"`
#'   return a different likelihood and are refused.
#' @param atol,rtol Absolute and relative solver tolerances.
#' @param on_error What to do about a solve that fails: `"penalize"`
#'   (the default) fills that group's rows with `penalty` and warns,
#'   naming the group; `"error"` stops instead, also naming it. Read
#'   "Failed solves" below first: on the automatic-differentiation tape
#'   most failures cannot be detected at all, so neither setting has the
#'   reach it appears to have during a fit.
#' @param penalty The filler value for `on_error = "penalize"`. It is on
#'   the scale of the nonlinear body's result, before the response
#'   link, so lower it for a model whose link exponentiates.
#' @param ... Further arguments for [deSolve::ode()].
#'
#' @return A numeric vector of length `nrow(data)` when one state is
#'   selected, otherwise a matrix with `nrow(data)` rows. On the
#'   automatic-differentiation tape both carry the `advector` class.
#'
#' @section Installation:
#' `frm_ode()` needs \pkg{RTMBode}, which is not on CRAN:
#' ```r
#' install.packages("RTMBode", repos = c(
#'   "https://kaskr.r-universe.dev",
#'   "https://cloud.r-project.org"))
#' ```
#'
#' @seealso [frm_ode_failures()] for the groups a penalty was written
#'   into, [bf()] for the nonlinear formula grammar, and
#'   `vignette("ode")` for a worked population pharmacokinetic model.
#'
#' @examples
#' # One-compartment oral pharmacokinetics with between-subject
#' # variability on the absorption and elimination rates.
#' #   dA/dt = -ka A            A(0) = dose
#' #   dC/dt =  ka A / V - ke C C(0) = 0
#' pk_dyn <- function(t, y, p) {
#'   "c" <- RTMB::ADoverload("c")
#'   list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
#' }
#'
#' set.seed(2026)
#' tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
#' n_id <- 6
#' dd <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt))),
#'                  time = rep(tt, n_id), dose = 100)
#' ka <- exp(rnorm(n_id, 0, 0.3))[as.integer(dd$id)]
#' ke <- exp(rnorm(n_id, log(0.2), 0.25))[as.integer(dd$id)]
#' dd$conc <- 100 * ka / (10 * (ka - ke)) *
#'   (exp(-ke * dd$time) - exp(-ka * dd$time)) + rnorm(nrow(dd), 0, 0.3)
#'
#' if (requireNamespace("RTMBode", quietly = TRUE)) {
#'   \donttest{
#'   fit <- frm(
#'     bf(conc ~ frm_ode(pk_dyn,
#'                       init   = list(dose, 0),
#'                       times  = time,
#'                       parms  = list(exp(lka), exp(lke), exp(lV)),
#'                       group  = id,
#'                       states = c("depot", "central"),
#'                       output = "central"),
#'        lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
#'       gaussian(),
#'     data = dd, start = list(beta = c(0, log(0.25), log(8))))
#'   fixef(fit)
#'   }
#' }
#' @export
frm_ode <- function(dynamics, init, times, parms, group = NULL,
                    output = NULL, states = NULL, t0 = 0,
                    method = "lsoda", atol = 1e-8, rtol = 1e-8,
                    on_error = c("penalize", "error"),
                    penalty = 1e6, ...) {
  ode_require_backend()
  # lexically scoped, so they must be established in the frame that
  # subassigns and concatenates - not in the caller (RTMB gotcha)
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")

  on_error <- match.arg(on_error)
  if (!is.function(dynamics)) {
    stop("`dynamics` must be a function(t, y, parms) returning the ",
         "derivatives", call. = FALSE)
  }
  if (length(method) != 1L || !is.character(method) ||
        !method %in% ode_adaptive_methods) {
    stop("`method` must name an adaptive integrator, one of: ",
         paste(ode_adaptive_methods, collapse = ", "),
         ". Fixed-step integrators (rk4, euler, ...) do not solve the ",
         "system to the tolerance the likelihood is defined at, so ",
         "they return a different likelihood, not a noisier one",
         call. = FALSE)
  }

  # The solve grid and the group membership decide the *structure* of
  # the tape, so they have to be known before it is built. An estimated
  # observation time or lag time is a different feature (and a harder
  # one), not an accident to let through as a cryptic coercion error.
  for (nm in c("times", "group", "t0")) {
    if (inherits(get(nm), "advector")) {
      stop("`", nm, "` is an estimated quantity. frm_ode() needs it as ",
           "data: it fixes the solve grid, which is built before the ",
           "tape. Estimated event times are not supported",
           call. = FALSE)
    }
  }
  times <- as.numeric(times)
  n_obs <- length(times)
  if (!n_obs) stop("`times` is empty", call. = FALSE)
  if (anyNA(times)) {
    stop("`times` contains NA", call. = FALSE)
  }

  if (is.null(group)) {
    gi <- rep(1L, n_obs)
    glab <- "1"
  } else {
    if (length(group) != n_obs) {
      stop("`group` has length ", length(group), " but `times` has ",
           n_obs, call. = FALSE)
    }
    gf <- if (is.factor(group)) droplevels(group) else factor(group)
    gi <- as.integer(gf)
    glab <- levels(gf)
  }
  groups <- split(seq_len(n_obs), gi)
  labels <- glab[as.integer(names(groups))]

  init_cols <- ode_columns(init, n_obs, "init")
  parm_cols <- ode_columns(parms, n_obs, "parms")
  t0_cols <- ode_columns(t0, n_obs, "t0")
  if (length(t0_cols) != 1L) {
    stop("`t0` must be a single column", call. = FALSE)
  }
  ode_check_constant(init_cols, groups, "init", labels)
  ode_check_constant(parm_cols, groups, "parms", labels)
  ode_check_constant(t0_cols, groups, "t0", labels)

  n_state <- length(init_cols)
  if (!is.null(states)) {
    if (length(states) != n_state) {
      stop("`states` names ", length(states), " states but `init` has ",
           n_state, " column", if (n_state == 1L) "" else "s",
           call. = FALSE)
    }
    states <- as.character(states)
  }
  if (n_state >= ode_state_warn) {
    ode_warn_once(
      paste0("states", n_state),
      "frm_ode() is solving a system of ", n_state, " states. Above ",
      "about ", ode_state_warn, " states the second-order derivative ",
      "path the Laplace approximation needs can return NaN gradients ",
      "(the failure is integrator-dependent and silent). Check the ",
      "gradient with diagnose() before trusting a fit this large")
  }

  # which states to return
  if (is.null(output)) {
    out_idx <- seq_len(n_state)
  } else if (is.character(output)) {
    if (is.null(states)) {
      stop("`output` is character, so `states` must name the states",
           call. = FALSE)
    }
    out_idx <- match(output, states)
    if (anyNA(out_idx)) {
      stop("`output` names a state that is not in `states`: ",
           paste(output[is.na(out_idx)], collapse = ", "),
           " (states: ", paste(states, collapse = ", "), ")",
           call. = FALSE)
    }
  } else {
    out_idx <- as.integer(output)
    if (anyNA(out_idx) || any(out_idx < 1L) || any(out_idx > n_state)) {
      stop("`output` must index states 1 to ", n_state, call. = FALSE)
    }
  }

  # deSolve needs list(dydt); a bare derivative vector is friendlier to
  # write and costs one wrapper
  func <- function(t, y, p) {
    r <- dynamics(t, y, p)
    if (is.list(r)) r else list(r)
  }

  col_at <- function(cl, i) if (length(cl) == 1L) cl else cl[i]

  outs <- lapply(out_idx, function(k) numeric(n_obs))
  failed <- character(0)

  for (g in seq_along(groups)) {
    idx <- groups[[g]]
    idx <- idx[order(times[idx])]
    i1 <- idx[1L]
    tstart <- as.numeric(col_at(t0_cols[[1L]], i1))
    if (times[idx[1L]] < tstart) {
      stop("group '", labels[[g]], "' has an observation time (",
           format(times[idx[1L]]), ") before t0 (", format(tstart),
           "); frm_ode() integrates forward from t0", call. = FALSE)
    }
    # deSolve requires times[1] to be the initial time; the duplicate
    # that an observation at t0 creates is tolerated, and the extra row
    # is dropped from the solution
    grid <- c(tstart, times[idx])
    y0 <- do.call(c, lapply(init_cols, col_at, i = i1))
    pv <- do.call(c, lapply(parm_cols, col_at, i = i1))

    gave_up <- NULL
    sol <- tryCatch(
      {
        # deSolve reports "integration was not successful" as a warning
        # and still returns a full-length matrix of whatever it reached,
        # so the warning is the only evidence that the numbers are not a
        # solution. It is left to propagate as well: during a fit it is
        # the user's signal, and there it is all we have.
        s <- withCallingHandlers(
          RTMBode::ode(y = y0, times = grid, func = func, parms = pv,
                       method = method, atol = atol, rtol = rtol, ...),
          warning = function(w) {
            if (grepl(ode_giveup_pattern, conditionMessage(w))) {
              gave_up <<- conditionMessage(w)
            }
          }
        )
        # a solver that gives up can also return a short matrix; that
        # must be a failure, not a length error later
        if (nrow(s) != length(grid)) {
          stop("the integrator returned ", nrow(s), " of ",
               length(grid), " requested time points", call. = FALSE)
        }
        # The remaining two checks work on numbers only: RTMB refuses
        # comparison on AD types, so on the tape a diverging trajectory
        # cannot be seen at all. See the "Failed solves" section of the
        # help page for what that means for on_error.
        if (!inherits(s, "advector")) {
          if (!is.null(gave_up)) {
            stop("the integrator did not converge: ", gave_up,
                 call. = FALSE)
          }
          if (!all(is.finite(s))) {
            stop("the integrator returned non-finite values",
                 call. = FALSE)
          }
        }
        s
      },
      error = function(e) e
    )

    if (inherits(sol, "condition")) {
      if (on_error == "error") {
        stop("frm_ode() failed to solve group '", labels[[g]], "': ",
             conditionMessage(sol), call. = FALSE)
      }
      failed <- c(failed, labels[[g]])
      for (k in seq_along(out_idx)) {
        z <- outs[[k]]
        z[idx] <- penalty
        outs[[k]] <- z
      }
      next
    }

    for (k in seq_along(out_idx)) {
      z <- outs[[k]]
      z[idx] <- sol[-1L, out_idx[k] + 1L]
      outs[[k]] <- z
    }
  }

  # A penalty is a value the caller did not ask for standing in for a
  # solution, so it is never silent. The log is set on every call, so a
  # later clean call clears a stale record.
  ode_failure_log$groups <- failed
  ode_failure_log$n_groups <- length(groups)
  ode_failure_log$penalty <- penalty
  ode_failure_log$when <- Sys.time()
  if (length(failed)) {
    warning("frm_ode(): the solve failed for ", length(failed),
            " of ", length(groups), " group",
            if (length(groups) == 1L) "" else "s", " (",
            paste(utils::head(failed, 10L), collapse = ", "),
            if (length(failed) > 10L) ", ..." else "",
            "). Their rows hold penalty = ", format(penalty),
            ", not a solution. If this came from predict(), simulate() ",
            "or residuals(), those rows of the result are the penalty ",
            "value. See the 'Failed solves' section of ?frm_ode for ",
            "what can and cannot be caught while fitting; ",
            "frm_ode_failures() repeats this.", call. = FALSE)
  }

  if (length(outs) == 1L) outs[[1L]] else do.call(cbind, outs)
}

#' Groups whose ODE solve failed in the last `frm_ode()` call
#'
#' Reports the groups that [frm_ode()] could not solve on its most
#' recent call, and whose rows therefore hold the `penalty` value rather
#' than a solution. `frm_ode()` also warns when that happens; this
#' function is for reading the record afterwards, for example after a
#' `predict()` whose warnings were suppressed.
#'
#' The record covers the last call only, from anywhere: a fit, a
#' `predict()`, or a direct call. It is reset by every `frm_ode()` call,
#' so a clean call clears it. See the "Failed solves" section of
#' [frm_ode()] for why a failure during fitting usually cannot be
#' recorded at all.
#'
#' @return `NULL` when the last call solved every group. Otherwise a
#'   list with `groups` (the labels that failed), `n_groups` (how many
#'   groups the call had), `penalty` (the value written into their
#'   rows), and `when`.
#' @seealso [frm_ode()]
#' @examples
#' # NULL until a solve fails
#' frm_ode_failures()
#' @export
frm_ode_failures <- function() {
  if (!length(ode_failure_log$groups %||% character(0))) return(NULL)
  list(groups = ode_failure_log$groups,
       n_groups = ode_failure_log$n_groups,
       penalty = ode_failure_log$penalty,
       when = ode_failure_log$when)
}

# ---------------------------------------------------------------------
# Frame-time structural guard.
#
# The within-group constancy of a *linear predictor* cannot be checked
# on the tape: RTMB refuses comparison on AD types, so frm_ode() cannot
# see the values it is handed. It can be checked exactly at frame
# assembly, where the design matrices and the grouping column are both
# in hand, and where the answer does not depend on the current parameter
# values. That matters: the failure mode this catches (a time-varying
# covariate on a dynamics parameter) is invisible to a value-based check
# at a start value of zero, which is the usual start value.

#' Every `frm_ode()` call inside an expression, matched to its formals.
#'
#' @noRd
find_ode_calls <- function(expr, out = list()) {
  if (is.call(expr)) {
    fn <- expr[[1L]]
    nm <- if (is.name(fn)) as.character(fn) else
      if (is.call(fn) && identical(as.character(fn[[1L]]), "::")) {
        as.character(fn[[3L]])
      } else ""
    if (identical(nm, "frm_ode")) {
      out <- c(out, list(
        tryCatch(match.call(frm_ode, expr), error = function(e) NULL)
      ))
    }
    for (i in seq_along(expr)) {
      if (i == 1L) next
      if (!missing_arg(expr[[i]])) out <- find_ode_calls(expr[[i]], out)
    }
  }
  out
}

missing_arg <- function(x) is.name(x) && !nzchar(as.character(x))

#' Refuse a dynamics parameter that varies inside a solve group.
#'
#' Walks the nonlinear body of every linear predictor for `frm_ode()`
#' calls, and for each nonlinear parameter appearing in that call's
#' `init` or `parms` argument checks that the parameter's design is
#' constant within the levels of the call's `group` column. Only a bare
#' symbol `group =` can be resolved against the model frame; anything
#' computed is left to the runtime check in [frm_ode()].
#'
#' @noRd
check_ode_constancy <- function(spec, linpreds, mf) {
  for (resp in spec$responses) {
    nl <- resp$dpars$mu$nl_body
    if (is.null(nl)) next
    calls <- Filter(Negate(is.null), find_ode_calls(nl))
    if (!length(calls)) next
    nlpars <- resp$nlpars %||% character(0)
    for (cl in calls) {
      gexpr <- cl[["group"]]
      if (is.null(gexpr) || !is.name(gexpr)) next
      gname <- as.character(gexpr)
      gv <- mf[[gname]]
      if (is.null(gv) || is.matrix(gv)) next
      gi <- as.integer(factor(gv))
      for (arg in c("init", "parms", "t0")) {
        aexpr <- cl[[arg]]
        if (is.null(aexpr)) next
        for (np in intersect(all.vars(aexpr), nlpars)) {
          lp <- linpreds[[linpred_key(resp$resp_name, np)]]
          if (is.null(lp)) next
          bad <- c(ode_varying_cols(lp$X, gi),
                   ode_varying_cols(lp$Z, gi))
          if (length(bad)) {
            stop("Nonlinear parameter '", np, "' is a dynamics input of ",
                 "frm_ode() but is not constant within '", gname,
                 "': ", paste(unique(bad), collapse = ", "),
                 ". frm_ode() solves one system per group and reads ",
                 "each dynamics input off the group's first row, so a ",
                 "term that varies within a group cannot enter the ",
                 "likelihood - the fit would leave its coefficient at ",
                 "the start value with an indefinite Hessian. Move the ",
                 "term to a nonlinear parameter that is not a dynamics ",
                 "input, or aggregate it to one value per group",
                 call. = FALSE)
          }
        }
      }
    }
  }
  invisible(TRUE)
}

#' Names of the columns of a design that vary inside some group.
#'
#' @noRd
ode_varying_cols <- function(M, gi) {
  if (is.null(M) || !nrow(M) || !ncol(M)) return(character(0))
  M <- as.matrix(M)
  cn <- colnames(M) %||% paste0("column ", seq_len(ncol(M)))
  # every row against its own group's first row: one comparison, and it
  # handles a sparse Z without a per-column scan over the whole b vector
  first <- stats::ave(seq_len(nrow(M)), gi, FUN = function(i) i[1L])
  d <- abs(M - M[first, , drop = FALSE])
  tol <- 1e-8 * max(1, max(abs(M)))
  cn[apply(d > tol, 2L, any)]
}
