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

#' Default for a `NULL` left side.
#'
#' frmtmb has one of these, but it is internal there and base R gained
#' `%||%` only in 4.4, later than this package's floor.
#'
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

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

# The two constants of the steady-state tail correction, both settled
# by the sweep in `dev/nss/nss-06-rule.R`.
#
# `ode_ss_rcap` is where the correction begins standing down, and it is
# applied in full below it. It is 0.9875 rather than 0.99 so that the
# LARGEST factor the smooth shape can apply, 99.02 at r = 0.99119, is
# the same bound a hard cap at 0.99 gave, 98.99: the shape is wider
# than the ramp it replaced, and left at 0.99 it would have raised that
# bound to 124. `ode_ss_noise` scales the integrator's own tolerance
# into the damping term, so a difference within a few multiples of the
# solver's noise sets no ratio. At 4 the sweep's worst row is 0.57
# against truncation's 0.61 and its median is 0.065 of it; at a hard
# cap and no damping the worst row is 36.
ode_ss_rcap <- 0.9875
ode_ss_noise <- 4
# Below this ratio the tail is at most 0.053 of the last difference, so
# gating it away is invisible and it lets the bottom end be closed
# smoothly instead of clamped.
ode_ss_rlow <- 0.05
# How much of the geometric tail the correction has to have summed
# before the numeric report is willing to quote that tail as a
# distance. Below it the stand-down has already said the ratio is not
# to be trusted, and the report falls back to the movement.
ode_ss_trust <- 0.5

# One-shot warning bookkeeping: the nl body is taped, so the R code
# below runs once or twice per fit, but re-fitting in the same session
# should not repeat a structural warning the user has already read.
ode_warned <- new.env(parent = emptyenv())

# The groups whose solve failed on the most recent frm_ode() call, so
# that a penalty written into a prediction can be read back after the
# fact. Reset by every call; see frm_ode_failures().
ode_failure_log <- new.env(parent = emptyenv())

#' Raise a solver warning at most once per key.
#'
#' A solve runs once per objective evaluation, so an unguarded warning
#' would repeat thousands of times in one fit. Keyed on `key`: the first
#' call warns and later calls are silent.
#'
#' @noRd
ode_warn_once <- function(key, ...) {
  if (!is.null(ode_warned[[key]])) return(invisible(NULL))
  ode_warned[[key]] <- TRUE
  warning(..., call. = FALSE)
}

#' Require the optional solver backend.
#'
#' @noRd
ode_has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

#' Refuse a `frm_ode()` call when no solver backend is installed.
#'
#' `RTMBode` and `deSolve` are both suggested, and `frm_ode()` needs
#' both. Returns invisibly when they are there, and otherwise errors,
#' naming the missing ones and where to get them.
#'
#' @noRd
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
ode_check_constant <- function(cols, groups, arg, labels,
                               who = "frm_ode()") {
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
             format(rng[2L]), "). ", who, " evaluates one system per ",
             "group and reads each dynamics input off the group's ",
             "first row, so a within-group covariate cannot enter the ",
             "likelihood.", call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# Dosing events.
#
# deSolve has an `events` argument, and RTMBode forwards `...` to it, but
# that route is unusable here for two separate reasons, both established
# in dev/ode-feasibility.md section 9:
#
#   1. It errors on the automatic-differentiation path. RTMBode hands
#      deSolve an unnamed state vector, deSolve bounds the event `var`
#      index by length(names(y)), and that bound is zero: "too many state
#      variables in 'event'; should be < 0". Upstream defect.
#   2. Worse, if that were fixed it would be silently wrong for two of
#      the three methods. RTMBode solves an AUGMENTED system whose tail
#      carries d(state)/d(parameter); an event row jumps the state and
#      leaves the sensitivity block alone. That is right for "add"
#      (adding a constant does not change a derivative) and wrong for
#      "replace" and "multiply", which measured 42% and 59% relative
#      gradient error against finite differences.
#
# So frm_ode() does not use deSolve events at all. It splits the
# integration at the event times and chains one solve per interval,
# carrying the end state forward and applying the jump in ordinary RTMB
# arithmetic. RTMBode is already differentiable through the initial
# state, so the adjoint is correct by construction for every method, and
# the dose amount is allowed to be an estimated quantity.
#
# "reset" is not a deSolve method. It is NONMEM's EVID = 3 and rxode2's
# evid = 3: every compartment is set to `value` (zero, in the records
# those two write), and the run continues from there. A reset row plus a
# dose row at the same instant is EVID = 4.
ode_event_methods <- c("add", "replace", "multiply", "reset")

#' Resolve a state selector (`output`, or an event's `state`) to positions.
#'
#' @noRd
ode_state_index <- function(x, states, n_state, what) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    if (is.null(states)) {
      stop("`", what, "` is character, so `states` must name the states",
           call. = FALSE)
    }
    i <- match(x, states)
    if (anyNA(i)) {
      stop("`", what, "` names a state that is not in `states`: ",
           paste(unique(x[is.na(i)]), collapse = ", "),
           " (states: ", paste(states, collapse = ", "), ")",
           call. = FALSE)
    }
    return(i)
  }
  i <- suppressWarnings(as.integer(x))
  if (anyNA(i) || any(i < 1L) || any(i > n_state)) {
    stop("`", what, "` must index states 1 to ", n_state, call. = FALSE)
  }
  i
}

#' Validate a dosing table and split it per group.
#'
#' Returns a list, one element per group label, holding that group's
#' events sorted by time, or `NULL` for a group with none. Everything
#' here is data: the table fixes the breakpoints of the solve, which are
#' chosen before the tape is built.
#'
#' @noRd
ode_split_events <- function(events, labels, n_state, states) {
  # A schedule that is not fixed at parse time - read from a file, or
  # built per fit - is passed as a function of no arguments and called
  # here. A data.frame held in a variable needs no such wrapper: it is
  # not a possible column, so a nonlinear body reads it from the formula
  # environment (drop_nl_lexical_datavars(), R/frame.R).
  if (is.function(events)) {
    events <- tryCatch(events(), error = function(e) {
      stop("`events` is a function and calling it failed: ",
           conditionMessage(e), call. = FALSE)
    })
    if (!is.data.frame(events)) {
      stop("`events` is a function, so it must return a data.frame of ",
           "doses; it returned ", class(events)[1L], call. = FALSE)
    }
  }
  if (inherits(events, "advector")) {
    stop("`events` is an estimated quantity. The event times fix where ",
         "the solve is split, which is decided before the tape is built, ",
         "so the table must be data. An estimated dose AMOUNT is ",
         "supported through `event_scale`", call. = FALSE)
  }
  if (!is.data.frame(events)) {
    stop("`events` must be a data.frame with columns time, value and ",
         "state (plus optional group, method and duration), one row per ",
         "dose", call. = FALSE)
  }
  if (!nrow(events)) {
    stop("`events` has no rows; pass NULL for a model without doses",
         call. = FALSE)
  }
  nms <- names(events)
  miss <- setdiff(c("time", "value"), nms)
  if (length(miss)) {
    stop("`events` is missing the ", paste(miss, collapse = " and "),
         " column", if (length(miss) > 1L) "s" else "",
         ". The columns are: group, time, state, value, method, ",
         "duration", call. = FALSE)
  }
  known <- c("group", "time", "state", "value", "method", "duration",
             "ii", "addl", "ss")
  extra <- setdiff(nms, known)
  if (length(extra)) {
    stop("`events` has unknown column", if (length(extra) > 1L) "s" else "",
         ": ", paste(extra, collapse = ", "),
         ". The columns are: ", paste(known, collapse = ", "),
         ". frm_ode() does not read NONMEM records: an `evid`/`amt`/`cmt` ",
         "table has to be reshaped to these names first",
         call. = FALSE)
  }

  time <- events[["time"]]
  if (!is.numeric(time) || anyNA(time) || any(!is.finite(time))) {
    stop("`events$time` must be finite and numeric", call. = FALSE)
  }
  value <- events[["value"]]
  if (!is.numeric(value) || anyNA(value) || any(!is.finite(value))) {
    stop("`events$value` must be finite and numeric. A dose that ",
         "depends on an estimated parameter goes in `event_scale`, not ",
         "here", call. = FALSE)
  }

  method <- if ("method" %in% nms) {
    m <- events[["method"]]
    if (is.factor(m)) m <- as.character(m)
    if (!is.character(m)) {
      stop("`events$method` must be one of ",
           paste(ode_event_methods, collapse = ", "), call. = FALSE)
    }
    bad <- setdiff(unique(m), ode_event_methods)
    if (length(bad)) {
      stop("`events$method` has unknown method",
           if (length(bad) > 1L) "s" else "", ": ",
           paste(bad, collapse = ", "), ". It must be one of ",
           paste(ode_event_methods, collapse = ", "), call. = FALSE)
    }
    m
  } else rep("add", nrow(events))

  # A reset sets every compartment, so it names no compartment. The
  # column is still required for every other row of a multi-state system.
  state <- if ("state" %in% nms) {
    st <- events[["state"]]
    if (is.factor(st)) st <- as.character(st)
    if (anyNA(st)) {
      if (any(is.na(st) & method != "reset")) {
        stop("`events$state` is NA on a row whose method is not ",
             "\"reset\". Only a reset names no compartment, because it ",
             "sets every one of them", call. = FALSE)
      }
      st[is.na(st)] <- if (is.character(st)) states[[1L]] else 1L
    }
    st
  } else {
    if (n_state != 1L && any(method != "reset")) {
      stop("`events` has no `state` column and the system has ", n_state,
           " states, so there is no state to dose. Name the compartment ",
           "each row goes into", call. = FALSE)
    }
    rep(1L, nrow(events))
  }
  state <- ode_state_index(state, states, n_state, "events$state")

  duration <- if ("duration" %in% nms) {
    dur <- events[["duration"]]
    if (!is.numeric(dur)) {
      stop("`events$duration` must be numeric", call. = FALSE)
    }
    dur[is.na(dur)] <- 0
    if (any(!is.finite(dur)) || any(dur < 0)) {
      stop("`events$duration` must be finite and not negative. Use 0 or ",
           "NA for an instantaneous dose", call. = FALSE)
    }
    if (any(dur > 0 & method != "add")) {
      stop("`events$duration` is positive on a row whose method is not ",
           "\"add\". An infusion delivers `value` at a constant rate ",
           "into the state, which is an addition; \"replace\" and ",
           "\"multiply\" are instantaneous only", call. = FALSE)
    }
    dur
  } else rep(0, nrow(events))

  # NONMEM's compact repetition. `ii` is the interdose interval, `addl`
  # the number of ADDITIONAL doses after the row's own. Both are pure
  # data preprocessing: the rows are written out below and nothing
  # downstream knows they were ever abbreviated.
  ii <- if ("ii" %in% nms) {
    v <- events[["ii"]]
    if (!is.numeric(v)) {
      stop("`events$ii` must be numeric: it is the interdose interval",
           call. = FALSE)
    }
    v[is.na(v)] <- 0
    if (any(!is.finite(v)) || any(v < 0)) {
      stop("`events$ii` must be finite and not negative. Use 0 or NA on ",
           "a row that is neither repeated nor at steady state",
           call. = FALSE)
    }
    v
  } else rep(0, nrow(events))

  addl <- if ("addl" %in% nms) {
    v <- events[["addl"]]
    if (!is.numeric(v)) {
      stop("`events$addl` must be numeric: it counts the doses that ",
           "follow the row's own", call. = FALSE)
    }
    v[is.na(v)] <- 0
    if (any(!is.finite(v)) || any(v < 0) || any(v != trunc(v))) {
      stop("`events$addl` must be a whole number and not negative",
           call. = FALSE)
    }
    as.integer(v)
  } else rep(0L, nrow(events))

  ss <- if ("ss" %in% nms) {
    v <- events[["ss"]]
    if (is.numeric(v)) v <- v != 0
    if (!is.logical(v)) {
      stop("`events$ss` must be TRUE/FALSE (or 1/0): it marks a row ",
           "whose dosing cycle has already reached steady state",
           call. = FALSE)
    }
    v[is.na(v)] <- FALSE
    v
  } else rep(FALSE, nrow(events))

  if (any(addl > 0L & ii <= 0)) {
    stop("`events$addl` is positive on a row whose `ii` is not. The ",
         "additional doses land at time + ii, 2 * ii, ..., so the ",
         "interval has to be given", call. = FALSE)
  }
  if (any(ss & ii <= 0)) {
    stop("`events$ss` is TRUE on a row whose `ii` is not positive. A ",
         "steady state is a state under a dose repeated every `ii`",
         call. = FALSE)
  }
  if (any(ss & method != "add")) {
    stop("`events$ss` is TRUE on a row whose method is not \"add\". A ",
         "steady state is reached by repeating a dose; \"replace\", ",
         "\"multiply\" and \"reset\" are not doses", call. = FALSE)
  }
  if (any(ss & duration > ii)) {
    stop("`events$ss` is TRUE on a row whose `duration` is longer than ",
         "its `ii`, so the infusions would overlap. Shorten the ",
         "duration, or lengthen the interval", call. = FALSE)
  }

  grp <- if ("group" %in% nms) {
    g <- events[["group"]]
    g <- if (is.factor(g)) as.character(g) else as.character(g)
    bad <- setdiff(unique(g), labels)
    if (length(bad)) {
      stop("`events$group` names ", length(bad), " group",
           if (length(bad) > 1L) "s" else "", " that ",
           if (length(bad) > 1L) "are" else "is",
           " not in `group`: ",
           paste(utils::head(bad, 5L), collapse = ", "),
           if (length(bad) > 5L) ", ..." else "", call. = FALSE)
    }
    g
  } else {
    # no group column: the same schedule for every subject
    NULL
  }

  tab <- data.frame(time = as.numeric(time), state = as.integer(state),
                    value = as.numeric(value),
                    method = as.character(method),
                    duration = as.numeric(duration),
                    ii = as.numeric(ii), ss = as.logical(ss),
                    grp = if (is.null(grp)) NA_character_ else grp,
                    row = seq_len(nrow(events)),
                    stringsAsFactors = FALSE)

  if (any(addl > 0L)) {
    tab <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) {
      if (addl[i] <= 0L) return(tab[i, , drop = FALSE])
      k <- 0:addl[i]
      x <- tab[rep(i, length(k)), , drop = FALSE]
      x$time <- tab$time[i] + k * ii[i]
      # the steady state is reached once, at the first of the repeats
      x$ss <- c(tab$ss[i], rep(FALSE, addl[i]))
      x
    }))
    tab$row <- seq_len(nrow(tab))
  }

  per <- if (is.null(grp)) {
    tab$grp <- NULL
    stats::setNames(rep(list(tab), length(labels)), labels)
  } else {
    s <- split(tab, factor(tab$grp, levels = labels))
    stats::setNames(lapply(labels, function(l) {
      x <- s[[l]]
      if (!is.null(x)) x$grp <- NULL
      x
    }), labels)
  }
  lapply(per, function(x) {
    if (is.null(x) || !nrow(x)) return(NULL)
    x <- x[order(x$time, x$row), , drop = FALSE]
    if (sum(x$ss) > 1L) {
      stop("`events` marks more than one row `ss = TRUE` for one group. ",
           "A run-in starts from an empty system, so a second one would ",
           "discard the first; write the later doses out instead",
           call. = FALSE)
    }
    # two rows on the same state at the same instant compose only when
    # both are additions; "replace" and "multiply" would depend on the
    # order of the rows, which is not a model the table can express. A
    # "reset" is exempt: it is always applied first, so a reset beside a
    # dose is NONMEM's EVID = 4 and its order is fixed, not ambiguous.
    q <- x[x$method != "reset", , drop = FALSE]
    key <- paste(q$time, q$state)
    dup <- key %in% key[duplicated(key)]
    if (any(dup & q$method != "add")) {
      i <- which(dup & q$method != "add")[1L]
      stop("`events` has more than one row for state ", q$state[i],
           " at time ", format(q$time[i]),
           " and one of them is \"", q$method[i],
           "\". Their order would decide the result, so the table is ",
           "ambiguous. Only repeated \"add\" rows compose",
           call. = FALSE)
    }
    rownames(x) <- NULL
    x
  })
}

#' Piecewise-constant time-varying inputs, as blocks of one group's rows.
#'
#' Returns `bps`, the times at which some `tv` column changes (with `t0`
#' first), and `vals`, one entry per `tv` column holding that column's
#' value on each of those intervals. `[bps[i], bps[i + 1])` carries
#' `vals[[j]][i]`: last observation carried forward, the convention
#' rxode2 uses for covariates, so a row's value is in force from that
#' row's time until the next change.
#'
#' The change points have to be data. When a `tv` column carries
#' estimated values there is nothing to compare - RTMB refuses comparison
#' on AD types - so `tv_break` names the data column whose changes mark
#' them.
#'
#' @noRd
ode_tv_blocks <- function(tv_cols, brk, idx, times, tstart, label) {
  n <- length(idx)
  tt <- times[idx]
  chg <- rep(FALSE, n)
  if (!is.null(brk)) {
    k <- as.character(brk[idx])
    k[is.na(k)] <- "\r<NA>"
    if (n > 1L) chg <- c(FALSE, k[-1L] != k[-n])
  } else {
    for (v in tv_cols) {
      if (length(v) == 1L) next
      x <- as.numeric(v[idx])
      if (n > 1L) {
        tol <- 1e-8 * max(1, max(abs(x)))
        chg <- chg | c(FALSE, abs(diff(x)) > tol)
      }
    }
  }
  if (n > 1L && any(chg & c(FALSE, diff(tt) == 0))) {
    j <- which(chg & c(FALSE, diff(tt) == 0))[1L]
    stop("group '", label, "' has two observations at time ",
         format(tt[j]), " that disagree about a `tv` value. A ",
         "time-varying input is a step function of time, so one time ",
         "cannot carry two values", call. = FALSE)
  }
  starts <- c(1L, which(chg))
  # the first block's value reaches back to t0; every later block starts
  # at the observation time where the value changed
  bps <- c(tstart, tt[which(chg)])
  vals <- lapply(tv_cols, function(v) {
    if (length(v) == 1L) rep(v, length(starts)) else v[idx[starts]]
  })
  # a plain-numeric tv column has to agree with the declared breaks,
  # otherwise the model solved is not the model the data describes
  if (!is.null(brk)) {
    first <- rep(FALSE, n)
    first[1L] <- TRUE
    blk <- cumsum(chg | first)
    for (j in seq_along(tv_cols)) {
      v <- tv_cols[[j]]
      if (inherits(v, "advector") || length(v) == 1L) next
      x <- as.numeric(v[idx])
      tol <- 1e-8 * max(1, max(abs(x)))
      if (any(abs(x - x[starts][blk]) > tol)) {
        stop("`tv` column ", j, " changes inside a block of group '",
             label, "' that `tv_break` says is constant. The break ",
             "column decides where the solve is split, so every `tv` ",
             "value has to hold over the whole block", call. = FALSE)
      }
    }
  }
  list(bps = bps, vals = vals)
}

#' Run a dosing cycle to (approximate) steady state.
#'
#' The convergence test a steady state really wants - repeat until the
#' cycle-start state stops moving - branches on a value, which the tape
#' cannot do: the number of solves has to be settled before the tape is
#' built. So the run-in is a FIXED number of cycles, `n_ss`, given as
#' data, and after `n` of them a linear system is still the accumulation
#' factor raised to `n` away from the limit.
#'
#' That tail is not unknown, though: it is a geometric series whose
#' ratio the run-in has already measured. Three successive cycle-start
#' states give `d1 = y[n-1] - y[n-2]` and `d2 = y[n] - y[n-1]`, and
#' their ratio is the slowest mode's per-cycle contraction, so the
#' remaining sum is `d2 * r / (1 - r)`. Adding it is ordinary arithmetic
#' on tape variables, so unlike the convergence test it RUNS DURING A
#' FIT, at the parameters the fit is currently at, and it costs no extra
#' solve. This is Aitken's extrapolation, and for linear kinetics it is
#' exact up to the second-slowest mode.
#'
#' Three things stop it from being applied blind, each measured in
#' `dev/nss-findings.md`:
#'
#' - a state with no steady state at all, an AUC compartment being the
#'   ordinary example, has a ratio of 1 and an unbounded tail. The
#'   ratio is read PER STATE so that such a state cannot corrupt the
#'   states that do settle;
#' - a system that has already converged has differences at the
#'   integrator's own noise level, and a ratio read off noise is noise.
#'   The denominator carries the noise floor, so a difference that small
#'   returns `r` near 0 and no correction;
#' - a component that is GROWING has a ratio above 1, and reading it as
#'   "nearly stopped" is the worst possible reading. The weight `w`
#'   falls to 0 as `r` reaches 1, so such a component is left exactly
#'   where truncation leaves it.
#'
#' Off the tape, where the values are ordinary numbers, the extrapolant
#' from the last three cycles is compared with the one from the three
#' before it, and the caller warns when they still disagree. That is a
#' distance to the LIMIT rather than the cycle-to-cycle movement the
#' warning used to report, which understated it by about
#' `1 / (lambda_z * ii)`.
#'
#' Returns the state at `t_ss` BEFORE the record's own dose - the trough
#' - because the caller applies that dose as an ordinary "add" row.
#' Every compartment starts at zero, which is NONMEM's and rxode2's
#' reading of a steady-state record.
#'
#' @noRd
ode_run_in <- function(run, pv, state, amt, ii, dur, t_ss, n_state,
                       dyn_rate, n_ss, ss_tol, ss_acc, label,
                       ss_extrapolate = TRUE, atol = 1e-8,
                       rtol = 1e-8) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  zero <- numeric(n_state)
  y <- zero
  step <- function(y, from, to, rate) {
    s <- run(y, c(from, to),
             if (is.null(dyn_rate)) pv else c(pv, rate), dyn_rate)
    s[nrow(s), 1L + seq_len(n_state)]
  }
  # the three cycle-start states before the last one, so the tail can be
  # measured from four iterates and its own extrapolation checked
  keep <- vector("list", 3L)
  for (k in seq_len(n_ss)) {
    j <- k - (n_ss - 3L)
    if (j >= 1L) keep[[j]] <- y
    a0 <- t_ss - (n_ss - k + 1L) * ii
    if (dur > 0) {
      rate <- zero
      rate[state] <- amt / dur
      y <- step(y, a0, a0 + dur, rate)
      if (dur < ii) y <- step(y, a0 + dur, a0 + ii, zero)
    } else {
      y[state] <- y[state] + amt
      y <- step(y, a0, a0 + ii, zero)
    }
  }
  ode_run_in_finish(keep, y, ss_tol, ss_acc, label, ss_extrapolate,
                    atol, rtol)
}

#' The geometric tail of a run-in, and how far it says the run-in is
#' from its limit.
#'
#' Split out of `ode_run_in()` so it can be tested on iterates a caller
#' supplies, including the ones that made the guards necessary: an
#' unbounded state, a system already inside the integrator's noise, and
#' a component that is growing rather than shrinking.
#'
#' @noRd
ode_ss_extrapolate <- function(ya, yb, yc, atol, rtol) {
  d1 <- yb - ya
  d2 <- yc - yb
  # A difference no larger than the integrator's own tolerance carries
  # no information about the ratio; damping the denominator by it sends
  # `r` to 0 there rather than to whatever the noise happens to say.
  del <- (ode_ss_noise * (atol + rtol * abs(yc)))^2
  r <- (d1 * d2) / (d1 * d1 + del)
  # pmin()/pmax() compare, which RTMB refuses on the tape; abs() is
  # taped, and min(a, b) = (a + b - |a - b|) / 2 is the same function
  lo <- function(x, a) (x + a + abs(x - a)) / 2
  hi <- function(x, a) (x + a - abs(x - a)) / 2
  # `u` is 1 while the contraction is convincing and 0 once the measured
  # ratio reaches 1, and `um` keeps `1 - r` out of a denominator: it is
  # `1 - r` while that is at least `h`, and `h` after, so the factor is
  # finite AT r = 1 instead of 0 times Inf.
  h <- 1 - ode_ss_rcap
  u <- (1 - r) / h
  uc <- lo(hi(u, 1), 0)
  um <- lo(u, 1)
  # `smooth` is the degree-7 smootherstep S: S(0) = 0, S(1) = 1, and its
  # first three derivatives vanish at both ends. `step` is S(x) / x,
  # which is the shape the STAND-DOWN needs: `step'(1)` is -1, and -1 is
  # exactly the slope that joins `r / (1 - r)` below the cap, which is
  # why that junction comes out C3 rather than merely continuous. A
  # cubic smoothstep would leave a kink at r = 1.
  smooth <- function(x) x^4 * (35 - x * (84 - x * (70 - 20 * x)))
  step <- function(x) x^3 * (35 - x * (84 - x * (70 - 20 * x)))
  # The gate closes the bottom end, where a ratio below 0 says the
  # differences alternate and there is no geometric tail to sum. It
  # spans [0, ode_ss_rlow], where the tail is at most rlow / (1 - rlow)
  # of the last difference. A gate has to SATURATE FLAT, so it is S and
  # not S(x) / x: with the latter the same -1 that is right for the
  # stand-down left a corner of 1 / (1 - rlow) at rlow, and made the
  # correction overshoot the tail by up to 1.2487x below it.
  gate <- smooth(lo(hi(r / ode_ss_rlow, 1), 0))
  fac <- (r / (um * h)) * step(uc) * gate
  # `f` and `d` are already on the tape; returning them lets the numeric
  # report say how much of the tail the stand-down declined to sum,
  # which is where the warning used to go quiet at large `n_ss`.
  list(y = yc + d2 * fac, r = r, f = fac, d = d2)
}

#' Apply the tail correction, and say how far the answer that is
#' RETURNED still is from the limit.
#'
#' `keep` holds the three cycle-start states before the last one and
#' `y` is the last, so four iterates are available whenever `n_ss` is at
#' least three, three at `n_ss` = 2 and two at `n_ss` = 1.
#'
#' The correction is built INSIDE the `ss_extrapolate` branch, not
#' before it. Built before it, `ss_extrapolate = FALSE` puts the whole
#' correction on the tape and throws the value away: 147 nodes, 69
#' percent of the outer tape at `n_ss` = 20, on an arm documented as
#' unchanged.
#'
#' @noRd
ode_run_in_finish <- function(keep, y, ss_tol, ss_acc, label,
                              ss_extrapolate, atol, rtol) {
  have <- !vapply(keep, is.null, TRUE)
  last <- y
  three <- have[2L] && have[3L]
  if (isTRUE(ss_extrapolate) && three) {
    y <- ode_ss_extrapolate(keep[[2L]], keep[[3L]], last, atol,
                            rtol)[["y"]]
  }
  # The convergence report is a comparison, so it exists on the numeric
  # path only. See the "Steady-state dosing" section of ?frm_ode for why
  # a fit cannot repeat it. Because there is no tape here, the report
  # may form the correction again without costing the FALSE arm a node.
  if (is.null(ss_acc) || !have[3L] || inherits(y, "advector") ||
        inherits(keep[[3L]], "advector")) {
    return(y)
  }
  ode_run_in_report(keep, last, y, have, three, ss_tol, ss_acc, label,
                    ss_extrapolate, atol, rtol)
}

#' How far the returned run-in state is from its limit, on the numeric
#' path.
#'
#' Every quantity here has to be one with NO HOLE IN IT, and that is the
#' whole difficulty. Each guard in `ode_ss_extrapolate()` works by
#' driving the correction to exactly zero, so a report built only out of
#' the correction is silent on precisely the states the guards exist
#' for. Measured when it was: three constructions, three silent, one of
#' them returning a value 0.2701 from its limit
#' (`dev/rev-nss/rev-nss-08-failopen.R`). `move`, the movement between
#' the last two cycle-start states, is the quantity the run-in has
#' always been able to compute and is what every branch falls back to.
#'
#' @noRd
ode_run_in_report <- function(keep, last, y, have, three, ss_tol,
                              ss_acc, label, ss_extrapolate, atol,
                              rtol) {
  scale <- max(1e-12, max(abs(as.numeric(y))))
  move <- abs(as.numeric(last) - as.numeric(keep[[3L]]))
  rel <- max(move) / scale
  stalled <- FALSE
  if (three) {
    ext <- ode_ss_extrapolate(keep[[2L]], keep[[3L]], last, atol, rtol)
    r <- as.numeric(ext[["r"]])
    dd <- abs(as.numeric(ext[["d"]]))
    # The WHOLE geometric tail, what the correction actually summed of
    # it, and what the stand-down declined. Where the ratio says there
    # is no tail the difference is zero and the `bad` floor below takes
    # over. The declined part is what made the warning fall silent at
    # large `n_ss`: with the correction partly stood down, successive
    # extrapolants agree with each other while both sit short of the
    # limit.
    fac <- as.numeric(ext[["f"]])
    whole <- r / (1 - r)
    # Where the correction has substantially STOOD DOWN the code has
    # already declared the ratio untrustworthy, so `whole` may not be
    # quoted either: for a state with no steady state it is 1e+10 times
    # the last difference, which is arithmetically true of a divergent
    # series and useless in a warning. There the model-free bound, the
    # movement itself, is what gets reported.
    modelled <- is.finite(whole) & r < 1 &
      abs(fac) >= ode_ss_trust * abs(whole)
    full <- ifelse(modelled, dd * abs(whole), move)
    undone <- ifelse(modelled, dd * abs(whole - fac), move)
    resid <- if (have[1L]) {
      prev <- ode_ss_extrapolate(keep[[1L]], keep[[2L]], keep[[3L]],
                                 atol, rtol)
      abs(as.numeric(ext[["y"]]) - as.numeric(prev[["y"]]))
    } else dd * abs(fac)
    # With only three iterates `resid` is the correction's own size, so
    # the truncated arm would count the same tail twice; the
    # extrapolated arm keeps it, where it is the only bound available on
    # the correction it just made.
    rel <- if (isTRUE(ss_extrapolate)) max(undone + resid) / scale
           else max(if (have[1L]) full + resid else full) / scale
    # A state whose movement is not already negligible and whose ratio
    # is not a contraction inside (0, 1) has no geometric model at all,
    # so neither `tail` nor `resid` describes it and its own movement is
    # the only bound left. A converged state fails the ratio test too -
    # an oral depot's differences are zero to the last bit, which reads
    # as a ratio of 0 - so the movement test comes first, or every oral
    # model would carry this floor.
    moving <- move > ss_tol * scale
    bad <- moving & !(r > 0 & r < 1)
    if (any(bad)) rel <- max(rel, max(move[bad]) / scale)
    stalled <- any(bad & r >= 1)
  }
  # A ratio at or above 1 says a state is not contracting between cycles
  # at all, so no `n_ss` reaches a steady state for it. That is its own
  # evidence and is reported on its own: inside the `rel` gate it was
  # suppressed by the same zero correction that suppressed `rel`.
  if (isTRUE(stalled) || (is.finite(rel) && rel > ss_tol)) {
    ss_acc$groups <- unique(c(ss_acc$groups, label))
    ss_acc$rel <- max(ss_acc$rel %||% 0,
                      if (is.finite(rel)) rel else 0)
    ss_acc$stalled <- isTRUE(ss_acc$stalled) || isTRUE(stalled)
  }
  y
}

#' Solve one group's system, split at its dosing events and its
#' piecewise-constant time-varying inputs.
#'
#' `run` performs one `RTMBode::ode()` call with the shared solver
#' settings and the failure checks. The state is carried across
#' breakpoints by hand, so the derivative of a later observation with
#' respect to a dose amount flows through the initial state of the
#' following segment, which RTMBode already differentiates exactly.
#'
#' A `tv` value rides the same machinery as an infusion rate: it is
#' appended to the parameter vector and held constant over one segment,
#' so it is a tape input and is differentiated exactly.
#'
#' @noRd
ode_solve_events <- function(run, y0, pv, tvals, ev, tstart, n_state,
                             dyn_rate, scale, tvb = NULL, n_ss = 0L,
                             ss_tol = 1e-6, ss_acc = NULL, label = "1",
                             ss_extrapolate = TRUE, atol = 1e-8,
                             rtol = 1e-8) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")

  if (is.null(ev)) {
    ev <- data.frame(time = numeric(0), state = integer(0),
                     value = numeric(0), method = character(0),
                     duration = numeric(0), ii = numeric(0),
                     ss = logical(0), row = integer(0),
                     stringsAsFactors = FALSE)
  }
  # scaling only means anything for an amount, so frm_ode() has already
  # refused a scale on a table carrying "replace" or "multiply"; a
  # "reset" is a state assignment rather than a dose, so it is exempt
  val <- ev$value * scale
  for (h in which(ev$method == "reset")) val[h] <- ev$value[h]
  ev_end <- ev$time + ev$duration
  inf <- ev$duration > 0
  tmax <- tvals[length(tvals)]

  # Every breakpoint is an exact double taken from `ev$time`,
  # `ev_end`, `tvb$bps` or `tvals`, so the equality tests below are
  # identity tests, not tolerance tests.
  bps <- sort(unique(c(tstart, ev$time[ev$time <= tmax],
                       ev_end[inf & ev_end <= tmax],
                       if (is.null(tvb)) NULL else tvb$bps[tvb$bps <= tmax],
                       tmax)))
  bps <- bps[bps >= tstart]
  n_seg <- length(bps) - 1L

  # last observation carried forward: the value in force on [a, b) is
  # the one the most recent change point set
  parms_at <- function(a) {
    if (is.null(tvb)) return(pv)
    ti <- findInterval(a, tvb$bps)
    # RTMB's c() refuses a NULL argument, and `parms` may legitimately
    # be empty when every input varies with time
    tvv <- lapply(tvb$vals, function(v) v[ti])
    do.call(c, if (is.null(pv)) tvv else c(list(pv), tvv))
  }

  # A steady-state record at t0 has to run before `init` is read: it
  # says what the state IS there, having been dosed this way for a long
  # time, so it replaces `init` rather than jumping from it. Doing it
  # here also covers a group whose only observation is at t0, where
  # there is no segment for the loop below to walk.
  ss_done <- rep(FALSE, nrow(ev))
  for (h in which(ev$ss & ev$time == tstart)) {
    y0 <- ode_run_in(run, parms_at(tstart), ev$state[h], val[h],
                     ev$ii[h], ev$duration[h], tstart, n_state, dyn_rate,
                     n_ss, ss_tol, ss_acc, label, ss_extrapolate,
                     atol, rtol)
    ss_done[h] <- TRUE
  }

  out <- lapply(seq_len(n_state), function(k) numeric(length(tvals)))
  # findInterval(left.open) sends t == tstart to 0: an observation at the
  # origin reads the initial state, before any event at the origin, the
  # same pre-dose reading deSolve gives at a dose time
  seg <- findInterval(tvals, bps, left.open = TRUE)
  at0 <- which(seg == 0L)
  if (length(at0)) {
    for (k in seq_len(n_state)) {
      z <- out[[k]]
      z[at0] <- y0[k]
      out[[k]] <- z
    }
  }

  y <- y0
  for (s in seq_len(n_seg)) {
    a <- bps[s]
    b <- bps[s + 1L]
    pvs <- parms_at(a)
    for (h in which(ev$ss & !ss_done & ev$time == a)) {
      y <- ode_run_in(run, pvs, ev$state[h], val[h], ev$ii[h],
                      ev$duration[h], a, n_state, dyn_rate, n_ss,
                      ss_tol, ss_acc, label, ss_extrapolate, atol, rtol)
      # An observation at the record's own time reads the run-in's
      # trough, overriding the value the preceding segment left, and it
      # is still pre-dose.
      at_a <- which(tvals == a)
      for (k in seq_len(n_state)) {
        z <- out[[k]]
        z[at_a] <- y[k]
        out[[k]] <- z
      }
    }
    # a reset comes first, so a reset beside a dose at one instant is
    # NONMEM's EVID = 4 rather than an ordering question
    for (h in which(ev$method == "reset" & ev$time == a)) {
      for (k in seq_len(n_state)) y[k] <- val[h]
    }
    for (h in which(ev$time == a & !inf & ev$method != "reset")) {
      k <- ev$state[h]
      y[k] <- switch(ev$method[h],
                     add = y[k] + val[h],
                     replace = val[h],
                     multiply = y[k] * val[h])
    }
    rate <- numeric(n_state)
    for (h in which(inf & ev$time <= a & ev_end > a)) {
      k <- ev$state[h]
      rate[k] <- rate[k] + val[h] / ev$duration[h]
    }
    rows <- which(seg == s)
    grid <- c(a, tvals[rows])
    if (grid[length(grid)] < b) grid <- c(grid, b)
    if (is.null(dyn_rate)) {
      sol <- run(y, grid, pvs, NULL)
    } else {
      sol <- run(y, grid, c(pvs, rate), dyn_rate)
    }
    if (length(rows)) {
      for (k in seq_len(n_state)) {
        z <- out[[k]]
        z[rows] <- sol[1L + seq_along(rows), k + 1L]
        out[[k]] <- z
      }
    }
    y <- sol[nrow(sol), 1L + seq_len(n_state)]
  }
  out
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
#' `tv` is the exception, and the next section is about it.
#'
#' # Time-varying inputs
#'
#' `tv` carries dynamics inputs that DO change inside a group, as long
#' as they change in steps: a creatinine clearance measured again at
#' each visit, a dose-dependent rate that switches at a protocol
#' amendment, a temperature held at one level and then another. The
#' solve is split at each change point, exactly as it is split at a dose
#' time, and each segment's dynamics see that segment's value.
#'
#' The values follow `parms` in the vector `dynamics` is handed, so with
#' `parms = list(ka, V)` and `tv = list(ke)` the derivative function
#' reads `p[1]`, `p[2]` and `p[3]`. Nothing in `dynamics` learns that
#' `p[3]` ever changed: within one segment it is an ordinary constant,
#' and because it is a parameter it is a tape input, so an estimated
#' time-varying value is differentiated exactly.
#'
#' The step function is **last observation carried forward**, which is
#' rxode2's convention for covariates. A row's value is in force from
#' that row's time until the next change, and the first row's value
#' reaches back to `t0`. A consequence worth knowing: the state is
#' continuous across a change point, because a covariate moves the
#' derivative and not the state, so an observation exactly at a change
#' point reads what the PRE-change dynamics produced. The difference
#' shows up only afterwards. (This is unlike a dose, which jumps the
#' state, and where the same observation reads the trough.)
#'
#' The values may be estimated; the change points may not. They decide
#' where the solve is split, which is settled before the tape is built.
#' When a `tv` column carries estimated values there is nothing to
#' compare - RTMB refuses comparison on AD types - so `tv_break` has to
#' name the data column whose changes within a group mark them:
#'
#' ```r
#' bf(conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
#'                   parms = list(exp(lka), exp(lV)),
#'                   tv = list(exp(lke)), tv_break = visit,
#'                   group = id, output = 2L),
#'    lka ~ 1, lV ~ 1, lke ~ 1 + crcl, nl = TRUE)
#' ```
#'
#' Here `lke` is a linear predictor that varies within a subject, which
#' `parms` would refuse. `tv_break = visit` says the value is constant
#' inside a visit, and a plain-numeric `tv` column that disagrees with
#' that is refused. When every `tv` column is plain data, `tv_break` can
#' be left out and the change points are read off the columns.
#'
#' A time-varying input other than a step function is not available. See
#' "Boundaries".
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
#' # Dosing events
#'
#' `events` is a data.frame of doses, one row per dose, with columns:
#'
#' - `group` (optional): which group the row applies to, matching the
#'   values of `group`. Leave the column out and the schedule applies to
#'   every group.
#' - `time`: when the dose happens. At or after `t0`.
#' - `state`: the state it goes into, by name (requires `states`) or by
#'   position, resolved exactly as `output` is. Optional for a
#'   one-state system, and ignored on a `"reset"` row, which names no
#'   compartment because it sets every one of them.
#' - `value`: how much.
#' - `method` (optional, default `"add"`): `"add"` puts `value` into the
#'   state, `"replace"` sets the state to `value`, `"multiply"` scales
#'   it - the \pkg{deSolve} event methods - and `"reset"` sets **every**
#'   state to `value`.
#' - `duration` (optional, default `0`): a positive value makes the row
#'   an infusion, delivering `value` at the constant rate
#'   `value / duration` over `[time, time + duration]`. Infusions must
#'   use `"add"`.
#' - `ii` (optional, default `0`): the interdose interval, for `addl`
#'   and `ss`.
#' - `addl` (optional, default `0`): how many further doses follow this
#'   one, at `time + ii`, `time + 2 * ii`, and so on. The rows are
#'   written out internally, so the table means exactly what the
#'   hand-expanded one means.
#' - `ss` (optional, default `FALSE`): the row's cycle has already
#'   reached steady state. See "Steady-state dosing".
#'
#' Inside a `bf(nl = TRUE)` body, name the table:
#'
#' ```r
#' doses <- data.frame(time = seq(12, 48, by = 12), state = "depot",
#'                     value = 100)
#' conc ~ frm_ode(pk_dyn, ..., events = doses)
#' ```
#'
#' A name in a nonlinear body is normally a request for a column of
#' `data`, and a column of that name still wins. A data.frame is not
#' something a column could hold, so `doses` is read from the formula
#' environment instead. Writing the table **inline**, or holding it in a
#' function of no arguments, works the same way and is what a schedule
#' read at fit time wants:
#'
#' ```r
#' conc ~ frm_ode(pk_dyn, ..., events = data.frame(
#'          time = seq(12, 48, by = 12), state = "depot", value = 100))
#'
#' schedule <- function() read.csv("doses.csv")
#' conc ~ frm_ode(pk_dyn, ..., events = schedule)
#' ```
#'
#' In NONMEM terms an `"add"` row is a dosing record (`evid = 1`) with
#' `amt = value` into `cmt = state`; a row with `duration` is the same
#' record with `rate = amt / duration`; `ii` and `addl` are spelled the
#' same way there and in rxode2. A `"reset"` row is `evid = 3`, which
#' both NONMEM and rxode2 read as "set every compartment to zero and
#' carry on" - write `value = 0` for that reading. `evid = 4`, a reset
#' followed by a dose, is a `"reset"` row and an `"add"` row at the same
#' time; the reset is always applied first, so the order of the two rows
#' in the table does not matter.
#'
#' `frm_ode()` does not read NONMEM column names, and there is no `evid`
#' column: observation rows are the rows of `data`, and dose rows are
#' the rows of `events`, which is a separate table. A NONMEM-shaped
#' dataset has to be split into the two.
#'
#' An observation at exactly a dose time reads the state **before** the
#' dose, which is the trough, matching both the \pkg{deSolve} convention
#' and the usual reading of a pre-dose sample. That includes an
#' observation at `t0` with a dose at `t0`: it reads `init`. The one
#' exception is an `ss` row, which is not a jump from anything.
#'
#' # Steady-state dosing
#'
#' `ss = TRUE` on an `events` row says the system has already been given
#' this dose every `ii` for long enough to settle. It is reached by
#' simulation, not by a closed form, so it works for any dynamics:
#' every compartment is set to zero (NONMEM's and rxode2's reading of a
#' steady-state record, and the reason an `ss` row overrides `init`),
#' and the cycle is then repeated `n_ss` times before the record's time.
#' The record's own dose is applied afterwards, so an observation at the
#' record's time reads the steady-state **trough**.
#'
#' ```r
#' # 100 into the depot every 12 hours, already at steady state at t = 0
#' data.frame(time = 0, state = "depot", value = 100, ii = 12, ss = TRUE)
#' ```
#'
#' A real steady state is the limit of repeating until the cycle-start
#' state stops moving, and that test branches on a value, which the tape
#' cannot do: the number of solves has to be fixed before the tape is
#' built. Truncating at `n_ss` cycles leaves a geometric tail. For
#' linear kinetics what is left after `n` cycles is the accumulation
#' factor to the power `n`, `exp(-n * lambda_z * ii)`, where `lambda_z`
#' is the SLOWEST disposition eigenvalue, so the shortfall is set by the
#' terminal half-life measured in dosing intervals and grows toward 1 as
#' that half-life grows. It moves estimates, not just trajectories: on
#' one 30-subject dataset simulated from the exact steady state,
#' truncating at 20 cycles rather than taking the limit moved `k21` by a
#' factor of 3.2 and `ke` by 11 percent.
#'
#' `ss_extrapolate = TRUE`, the default, **sums that tail instead of
#' dropping it**. The run-in already computes the last cycle-start
#' states, and their successive differences give the per-cycle
#' contraction `r`, so the sum that is left is `d * r / (1 - r)`. That
#' is arithmetic on tape variables, so it runs during a fit, at whatever
#' parameters the fit has reached, and it costs no extra solve. The
#' ratio is read one state at a time, damped by the integrator's own
#' tolerance, and stood down as `r` approaches 1. Each of the three is
#' needed; `dev/nss-findings.md` gives the construction that breaks the
#' rule without it.
#'
#' Measured worst over the run-in states against the exact limit, on a
#' two-compartment oral model at `n_ss = 20` and `atol = rtol = 1e-8`:
#'
#' | terminal half-life | `ii` | truncated | extrapolated |
#' | --- | --- | --- | --- |
#' | 23 h | 8 | 8.5e-03 | 7.6e-12 |
#' | 107 h | 24 | 4.5e-02 | 4.2e-12 |
#' | 107 h | 12 | 2.1e-01 | 7.8e-11 |
#' | 265 h | 24 | 2.8e-01 | 4.9e-11 |
#' | 670 h | 24 | 6.1e-01 | 5.0e-10 |
#'
#' The column above is worst over the run-in STATES; NEWS.md quotes the
#' same designs worst over one dosing interval on the observable, which
#' is the smaller number, 1.2e-02 rather than 4.5e-02 at 107 hours.
#' What is left after the correction is the integrator's own tolerance
#' **amplified by about `1 / (1 - r)`**, so it grows with the terminal
#' half-life: measured in units of `atol = rtol`, 0.3 to 1.4 on the five
#' rows above and 62.8 at a 1655 hour half-life dosed daily.
#'
#' # Where the correction is worse than truncating
#'
#' The ratio the run-in reads is a weighted mean of the cycle map's
#' modes. Where every weight has the same sign that mean lies between
#' the slowest and the fastest, so the correction can only undershoot.
#' Where two weights have opposite signs it lies OUTSIDE their range and
#' the correction overshoots, and opposite signs are the normal
#' arrangement in an oral model, which is why a concentration rises
#' before it falls. The hazard is therefore `ka` near `lambda_z`:
#' **flip-flop kinetics**, which extended-release and depot
#' formulations are written to produce.
#'
#' Measured over 338 two-compartment oral cycle maps, `lambda_z * ii`
#' from 0.02 to 3 and `ka / lambda_z` from 0.1 to 50: the correction
#' loses on **8** of them, worst by a factor of **2.27**, and every
#' losing case has `lambda_z * ii` = 0.05 with `ka / lambda_z` of 1.5 or
#' 2. It is bounded, and the bound is what makes the default defensible:
#' the smallest error truncation leaves on a losing case is **0.40**,
#' the largest it leaves on a winning one is 0.96, and the warning below
#' fires on every losing case in both arms. The correction never loses
#' where truncation was usable, and where it loses the user is told.
#' Through `frm_ode()` on the worst such model, a 333 hour terminal
#' half-life with `ka / lambda_z` = 1.5 dosed daily: 5.2e-01 truncated
#' against 8.0e-01 extrapolated, both warned.
#'
#' The correction is exact only up to the second-slowest mode for other
#' reasons too. On a Michaelis-Menten system deep in its saturated
#' regime it improved a 1.1e-01 shortfall to 7.2e-03 and no further, and
#' on a system with two modes of nearly the same rate it removes their
#' combination and leaves what is left of the other.
#'
#' # The stand-down, and the objective's smoothness
#'
#' The correction is applied in full while the measured ratio is between
#' 0.05 and 0.9875, gated to nothing below the first and stood down to
#' nothing as it reaches 1. The gate is the degree-7 smootherstep, whose
#' first three derivatives vanish at both ends; the stand-down is that
#' same shape divided by its argument, degree 6, whose slope at 1 is -1,
#' which is exactly the slope that joins `r / (1 - r)` below the cap.
#' Measured rather than asserted, **the objective has no corner in the
#' parameters** at any of the four junctions: a one-sided first
#' difference there reads the curvature `2 / (1 - r)^3` to five figures,
#' the same law it reads at ratios with no junction at all. That
#' matters because an optimizer walks into this region: a fit whose
#' `k21` runs to zero drives `lambda_z` to zero and `r` to 1, and it
#' crosses 0.9875 on the way. A hard cap there left the two one-sided
#' derivatives at -0.09 and +150, not converging as the bracket
#' tightened.
#'
#' Measured on the factor itself, as the gap between the two one-sided
#' difference quotients: at `r` = 0 it falls to exactly 0 by a bracket
#' of 1e-4; at `r` = 1 it falls like the square of the bracket; and at
#' `r` = 0.05 and `r` = 0.9875 the gap divided by the bracket is 2.333
#' and 1.024e+06, which are `2 / (1 - r)^3`, the factor's own second
#' derivative, to four figures. A one-sided first difference carries an
#' error of exactly that size, so at those two junctions the probe is
#' reading curvature and there is no discontinuity left in it to find.
#'
#' `ss_extrapolate = FALSE` restores the truncated run-in exactly, bit
#' for bit, and with the base commit's tape. Then the shortfall is
#' `exp(-n_ss * lambda_z * ii)` again, and **choose `n_ss` so that
#' `n_ss * lambda_z * ii` is at least 20**, which puts it at 2e-09. On
#' the 107 hour, `ii = 24` model that is `n_ss = 134`, not 20, and the
#' run-in is most of the solve count.
#'
#' A very long run-in has a second limit that raising `n_ss` cannot
#' pass. The cycles are chained solves, so the integrator's own error
#' accumulates across all `n_ss + 1` of them: at `atol = rtol = 1e-8`
#' and `n_ss` = 1000 it contributes about 2e-05, which is larger than
#' the run-in shortfall it was raised to remove. Past a few hundred
#' cycles, tighten `atol` and `rtol` as well or the extra cycles buy
#' nothing.
#'
#' # What frm_ode() will tell you
#'
#' Off the tape - a direct call, [predict()], [simulate()], a body
#' holding no estimated parameter - `frm_ode()` reports how far the
#' value IT RETURNED still is from the limit, and warns when that is
#' more than `ss_tol`. **That is a distance to the limit**, which the
#' warning before 0.4.0 was not: it reported the cycle-to-cycle
#' movement, which understates the distance by about
#' `1 / (lambda_z * ii)`, measured at 3.1x, 6.6x and 10.8x as the
#' half-life grows, so it understated it most exactly where the error
#' was largest. It also says so by name when a state is not contracting
#' between cycles at all, which means no `n_ss` settles it.
#'
#' **Treat it as a detector and not as a measurement.** It is built out
#' of the same geometric model the correction is, so where that model
#' is poor the number is poor with it, and it is bounded below by the
#' cycle-to-cycle movement rather than being an estimate in its own
#' right. On a state whose ratio is not a contraction it reports that
#' movement, which can be far from the distance to a limit that does
#' not exist. What it is good for is deciding whether to look, not how
#' much to trust the third digit.
#'
#' On 31 schedules classified against the exact limit there is no false
#' alarm and no miss in either arm, and the set now reaches both ends
#' of the range: three flip-flop models, `n_ss` of 1, 2 and 4, and a
#' 1653 hour half-life at `n_ss` of 650, 1000 and 2000, which is the
#' region the paragraph above sends a long-half-life user to.
#'
#' During a fit the check still cannot run at all. For a linear
#' compartment model [frm_lincmt()] sums the series in closed form and
#' has nothing to truncate.
#'
#' The cost is `n_ss` extra solves per group (two per cycle for an
#' infusion), so a steady-state population fit is several times a plain
#' one, and `ss_extrapolate` does not change that count. One `ss` row
#' per group is allowed; write later doses out with `ii` and `addl`.
#'
#' The doses are not handed to \pkg{deSolve} as events. `frm_ode()`
#' splits the integration at the event times and chains one solve per
#' interval, carrying the state across the break itself. That is a
#' correctness requirement, not a style choice: \pkg{RTMBode} solves an
#' augmented system carrying the derivatives of the states with respect
#' to the parameters, and a \pkg{deSolve} event jumps the state without
#' jumping those derivatives, which gives a wrong gradient for
#' `"replace"` and `"multiply"` (measured at 42% and 59% relative error).
#' Splitting the solve is exact for every method. It costs one solve per
#' dosing interval per group, so an intensively dosed design is
#' proportionally slower.
#'
#' # Doses that depend on a parameter
#'
#' `events` is data: `value` is a numeric column, so it cannot hold an
#' estimated quantity. `event_scale` is the way in. It is one value per
#' observation (constant within group, like `init` and `parms`), and it
#' multiplies the `value` of every event in that group, so a
#' bioavailability written as a nonlinear parameter estimates a dose
#' scale that carries covariates and random effects like any other:
#'
#' ```r
#' frm_ode(pk_dyn, init = list(0, 0), times = time, parms = ...,
#'         group = id, events = doses, event_scale = plogis(logitF))
#' ```
#'
#' Because scaling only makes sense for a dose, `event_scale` is refused
#' on a table containing `"replace"` or `"multiply"` rows. A `"reset"`
#' row is allowed beside scaled doses and is simply not scaled: it names
#' a level for the states, not an amount.
#'
#' # One solve, a different state per row
#'
#' `output` normally selects states. Given one value per row instead, it
#' says which state THAT row reads, and the whole group still takes one
#' solve. This is the spelling for a parent and its metabolite reported
#' in one assay column, or for two species of a predator-prey series
#' stacked long:
#'
#' ```r
#' bf(y ~ frm_ode(pm_dyn, init = list(dose, 0, 0), times = time,
#'                parms = list(exp(lka), exp(lkm), exp(lke)),
#'                group = id, output = cmt),          # cmt is 2 or 3
#'    lka ~ 1, lkm ~ 1, lke ~ 1, nl = TRUE)
#' ```
#'
#' The two readings are told apart by length: a per-row selection needs
#' one value per row AND more rows than the system has states. A
#' selection of states is shorter than that. The per-row values are
#' resolved by name or position exactly as a scalar `output` is, and
#' they have to be data - a state index cannot be differentiated.
#'
#' # Boundaries
#'
#' Time-varying input has to be piecewise constant, and the reason is
#' worth knowing. \pkg{RTMBode} tapes `dynamics` once, so `t` is an
#' automatic-differentiation value inside it. A branch on time
#' (`if (t < t_end) rate else 0`) raises "Comparison is generally unsafe
#' for AD types", and an `approxfun()` forcing table silently returns
#' the value at the taping point instead of failing. Smooth arithmetic in
#' `t` is fine. A step function belongs in `tv`, or in `events` as an
#' infusion, where it is carried as a parameter over each interval and
#' differentiated exactly; \pkg{deSolve}'s own `forcings` argument is not
#' reachable, because \pkg{RTMBode}'s compiled derivative shim has no
#' forcing hook. Linear interpolation between measured covariate values,
#' which rxode2 offers, is not available: it would need the value inside
#' the segment to depend on `t`.
#'
#' Estimated event times, lag times, inter-dose intervals and `tv`
#' change points are not supported: they decide where the solve is
#' split, which is settled before the tape is built.
#'
#' `predict(se.fit = TRUE)` is not available for a nonlinear predictor,
#' including one containing `frm_ode()`; request a nonlinear parameter
#' with `predict(dpar = )` instead.
#'
#' @param dynamics A function `function(t, y, parms)` giving the
#'   derivatives, following the \pkg{deSolve} convention: `t` is the
#'   scalar time, `y` the state vector, `parms` the parameter vector,
#'   and the return value is `list(dydt)`. A bare derivative vector is
#'   also accepted. Index `y` and `parms` by position. RTMB's
#'   tape-safe `c()`, `[<-` and `diag<-` are put in scope
#'   automatically, so the function needs no
#'   `"c" <- RTMB::ADoverload("c")` boilerplate; only a HELPER it
#'   calls, defined elsewhere, still needs its own (lexical scope does
#'   not travel into other functions).
#' @param init Initial states, one column per state: a list, a matrix,
#'   or a single vector for a one-state system. Each column is either
#'   one value per observation (constant within group) or one value
#'   shared by every group. The number of columns sets the number of
#'   states.
#' @param times Observation times, one per row of the data.
#' @param parms Dynamics parameters, one column per parameter, in the
#'   order `dynamics` expects them. Same shape rules as `init`. May be
#'   empty when every input is given in `tv`.
#' @param group Grouping vector naming the unit that owns one system,
#'   one value per row. `NULL` treats the whole data as one group.
#' @param output Which states to return: `NULL` (the default) returns
#'   every state, an integer or character vector selects some. A single
#'   selected state is returned as a vector, otherwise a matrix with one
#'   column per selected state. Character selection requires `states`.
#'   One value per row - which needs more rows than there are states -
#'   instead says which state each row reads, over one shared solve; see
#'   "One solve, a different state per row".
#' @param states Optional state names, one per column of `init`.
#' @param t0 Initial time, a scalar or one value per row (constant
#'   within group). Every observation time must be at or after it.
#' @param events Optional dosing table: a data.frame with columns
#'   `time`, `value` and `state`, and optional `group`, `method`,
#'   `duration`, `ii`, `addl` and `ss`, or a function of no arguments
#'   returning one. See "Dosing events" below. `NULL` (the default) is a
#'   model driven only by its initial conditions.
#' @param event_scale A multiplier on every `events$value`, one value per
#'   observation (constant within group) or a single value shared by
#'   every group. This is the one estimated quantity that can reach a
#'   dose, and it is how a bioavailability is written. Only for a table
#'   whose rows are all `method = "add"` or `"reset"`; a reset is not
#'   scaled.
#' @param tv Dynamics inputs that vary with time, one column per input,
#'   appended to `parms` in the vector `dynamics` receives. Unlike
#'   `parms` these may vary within a group, as a step function of time:
#'   the solve is split at each change point. See "Time-varying
#'   inputs".
#' @param tv_break A data column, one value per row, whose changes
#'   within a group are the change points of every `tv` column. Required
#'   when a `tv` column carries estimated values, because those cannot
#'   be compared. Optional otherwise, when the change points are read
#'   off the `tv` columns themselves.
#' @param n_ss How many dosing cycles a steady-state run-in simulates
#'   before an `events` row marked `ss = TRUE`. The steady state is
#'   approached, not solved for; see "Steady-state dosing" for the size
#'   of the approximation.
#' @param ss_tol Distance to the steady state, relative to the state's
#'   own scale, that `frm_ode()` will accept without warning. Checked on
#'   the numeric path only.
#' @param ss_extrapolate Sum the geometric tail that truncating the
#'   run-in at `n_ss` cycles leaves, using the contraction ratio the
#'   run-in itself measures. `TRUE`, the default, is exact for linear
#'   kinetics up to the second-slowest mode and costs no extra solve,
#'   though it does add about 170 nodes to the tape. `FALSE` truncates,
#'   which is what versions before 0.4.0 did, bit for bit. Read
#'   "Where the correction is worse than truncating" before assuming
#'   the default is always the better one. See "Steady-state dosing".
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
#'   selected, or when `output` selects one state per row, otherwise a
#'   matrix with `nrow(data)` rows. On the automatic-differentiation
#'   tape both carry the `advector` class.
#'
#' @section Installation:
#' `frm_ode()` needs \pkg{RTMBode}, which is not on CRAN:
#' ```r
#' install.packages("RTMBode", repos = c(
#'   "https://kaskr.r-universe.dev",
#'   "https://cloud.r-project.org"))
#' ```
#'
#' @section Sampling an ODE fit:
#' `frmtmb.sample::frm_sample()` and `frmtmb.sample::as_tmbstan()` both
#' refuse a model that contains `frm_ode()`, and the refusal is a
#' registered row of the compatibility registry
#' ([frmtmb::frm_compat()]) rather than a guess. Both doors are guarded,
#' and the second one was the worse of the two: it returned an empty
#' `stanfit` with no error at all, and the abort behind it left rstan
#' unusable for the rest of the R session.
#'
#' The defect is upstream and is reproduced in
#' `dev/upstream/rtmbode-issues.md` without frmtmb, in bare RTMB,
#' RTMBode and deSolve. `RTMBode::ode()` calls `deSolve::ode()` with no
#' guard. deSolve answers an extreme parameter either by raising
#' `illegal input detected before taking any integration steps` or by
#' returning fewer rows than were asked for, which the fixed-length
#' adjoint node reports as `Wrong output length`. An optimizer shortens
#' its step and a sampler rejects its proposal only when a failure
#' arrives as `NaN`; as an error it is fatal. Stan reaches those
#' parameters by construction, because its first warmup step is of size
#' 1 on the unconstrained scale, so the chain aborts at iteration 1 even
#' when it starts at the fitted optimum. The abort also crosses Stan's
#' C++ boundary and leaves rstan's nested autodiff arena unbalanced for
#' the rest of the R session.
#'
#' Two ways forward. Fit by maximum likelihood and read the Wald
#' intervals, which is what the rest of this page describes; or apply
#' the three-patch series in `dev/upstream/patches/` to RTMBode, after
#' which the same models sample. The refusal disappears on its own once
#' a fixed RTMBode is released and this package drops the row.
#'
#' @seealso [frm_ode_failures()] for the groups a penalty was written
#'   into, [frmtmb::bf()] for the nonlinear formula grammar, and
#'   `vignette("ode")` for a worked population pharmacokinetic model.
#'
#' @examples
#' # One-compartment oral pharmacokinetics with between-subject
#' # variability on the absorption and elimination rates.
#' #   dA/dt = -ka A            A(0) = dose
#' #   dC/dt =  ka A / V - ke C C(0) = 0
#' pk_dyn <- function(t, y, p) {
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
#'
#'   # Repeated dosing: 100 into the depot every 12 hours. The dose at
#'   # time 0 is the initial condition, the rest are events.
#'   doses <- data.frame(time = c(12, 24, 36), state = "depot",
#'                       value = 100)
#'   frm_ode(pk_dyn, init = list(100, 0), times = c(6, 18, 30, 42),
#'           parms = list(1, 0.2, 10), states = c("depot", "central"),
#'           output = "central", events = doses)
#'
#'   # The same schedule written compactly, and already at steady state
#'   frm_ode(pk_dyn, init = list(0, 0), times = c(6, 18, 30, 42),
#'           parms = list(1, 0.2, 10), states = c("depot", "central"),
#'           output = "central",
#'           events = data.frame(time = 0, state = "depot", value = 100,
#'                               ii = 12, ss = TRUE))
#'
#'   # An elimination rate that doubles after hour 6, carried forward
#'   # from the row it appears on. `tv` values follow `parms`, so this
#'   # dynamics reads ka at p[1], V at p[2] and the time-varying ke at
#'   # p[3].
#'   pk_tv <- function(t, y, p) {
#'     list(c(-p[1] * y[1], p[1] * y[1] / p[2] - p[3] * y[2]))
#'   }
#'   tt <- c(1, 3, 6, 9, 12)
#'   frm_ode(pk_tv, init = list(100, 0), times = tt, parms = list(1, 10),
#'           tv = list(ifelse(tt < 6, 0.2, 0.4)),
#'           states = c("depot", "central"), output = "central")
#' }
#' @export
frm_ode <- function(dynamics, init, times, parms = list(), group = NULL,
                    output = NULL, states = NULL, t0 = 0,
                    events = NULL, event_scale = 1,
                    tv = NULL, tv_break = NULL,
                    n_ss = 20L, ss_tol = 1e-6, ss_extrapolate = TRUE,
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
  # user code runs with the AD overloads in scope, so a dynamics
  # function need not carry the ADoverload boilerplate itself
  dynamics <- frmtmb::frmtmb_ad_overload(dynamics)
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
  for (nm in c("times", "group", "t0", "tv_break")) {
    if (inherits(get(nm), "advector")) {
      stop("`", nm, "` is an estimated quantity. frm_ode() needs it as ",
           "data: it fixes the solve grid, which is built before the ",
           "tape. Estimated event times are not supported",
           call. = FALSE)
    }
  }
  if (length(n_ss) != 1L || !is.numeric(n_ss) || is.na(n_ss) ||
        n_ss < 1 || n_ss != trunc(n_ss)) {
    stop("`n_ss` must be one whole number, at least 1: it is how many ",
         "dosing cycles a steady-state run-in simulates", call. = FALSE)
  }
  n_ss <- as.integer(n_ss)
  if (length(ss_extrapolate) != 1L || !is.logical(ss_extrapolate) ||
        is.na(ss_extrapolate)) {
    stop("`ss_extrapolate` must be TRUE or FALSE", call. = FALSE)
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
  # a system whose only inputs vary with time has no constant `parms`
  parm_cols <- if (is.null(parms) ||
                     (is.list(parms) && !inherits(parms, "advector") &&
                        !length(parms))) {
    if (is.null(tv)) {
      stop("`parms` is NULL or empty and `tv` was not given, so the ",
           "dynamics has no parameters at all. Pass the constants in ",
           "`parms`, the time-varying ones in `tv`", call. = FALSE)
    }
    list()
  } else ode_columns(parms, n_obs, "parms")
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

  # Which states to return. `output` normally selects states, but one
  # value per row selects a DIFFERENT state per row - a parent and its
  # metabolite measured on the same assay, a two-species series stacked
  # long - and then the solve is shared and the selection is a gather.
  # The two readings are told apart by length: a per-row selection needs
  # more rows than the system has states.
  out_sel <- NULL
  if (!is.null(output) && length(output) == n_obs && n_obs > n_state) {
    if (inherits(output, "advector")) {
      stop("`output` is an estimated quantity. It picks which state each ",
           "row reads, which is data: a state index cannot be ",
           "differentiated", call. = FALSE)
    }
    out_sel <- ode_state_index(output, states, n_state, "output")
    out_idx <- seq_len(n_state)
  } else {
    out_idx <- if (is.null(output)) seq_len(n_state) else
      ode_state_index(output, states, n_state, "output")
  }

  # Piecewise-constant time-varying inputs. Unlike `parms` these MAY
  # vary within a group: the solve is split at each change point, and
  # each segment's dynamics see that segment's value as an extra
  # parameter. The values may be estimated; the change points may not.
  tv_cols <- NULL
  n_tv <- 0L
  brk <- NULL
  if (!is.null(tv)) {
    tv_cols <- ode_columns(tv, n_obs, "tv")
    n_tv <- length(tv_cols)
    if (!is.null(tv_break)) {
      if (length(tv_break) != n_obs) {
        stop("`tv_break` has length ", length(tv_break), " but `times` ",
             "has ", n_obs, "; it is one value per observation",
             call. = FALSE)
      }
      brk <- tv_break
    } else if (any(vapply(tv_cols, function(v) inherits(v, "advector"),
                          TRUE))) {
      stop("`tv` carries an estimated quantity, so its change points ",
           "cannot be read off its values: RTMB refuses comparison on ",
           "AD types. Pass `tv_break`, the data column whose changes ",
           "within a group are the change points", call. = FALSE)
    }
  } else if (!is.null(tv_break)) {
    stop("`tv_break` was given but `tv` was not; there is nothing for ",
         "it to split", call. = FALSE)
  }

  # dosing events, and the one estimated quantity allowed to reach them
  ev_by_group <- NULL
  scale_cols <- NULL
  if (!is.null(events)) {
    ev_by_group <- ode_split_events(events, labels, n_state, states)
    has_infusion <- any(vapply(ev_by_group, function(x)
      !is.null(x) && any(x$duration > 0), TRUE))
    has_ss <- any(vapply(ev_by_group, function(x)
      !is.null(x) && any(x$ss), TRUE))
    if (!identical(event_scale, 1)) {
      meths <- unique(unlist(lapply(ev_by_group, function(x) x$method)))
      # a "reset" is a state assignment, not an amount, so it simply is
      # not scaled; "replace" and "multiply" would change meaning
      if (!all(meths %in% c("add", "reset"))) {
        stop("`event_scale` scales a dose, so it applies only to ",
             "\"add\" rows, but `events` also has ",
             paste(setdiff(meths, c("add", "reset")), collapse = " and "),
             " rows. Scaling those would change what they mean. Split ",
             "the model, or fold the scale into the state itself",
             call. = FALSE)
      }
      scale_cols <- ode_columns(event_scale, n_obs, "event_scale")
      if (length(scale_cols) != 1L) {
        stop("`event_scale` must be a single column", call. = FALSE)
      }
      ode_check_constant(scale_cols, groups, "event_scale", labels)
    }
  } else {
    has_infusion <- FALSE
    has_ss <- FALSE
    if (!identical(event_scale, 1)) {
      stop("`event_scale` was given but `events` was not; there is ",
           "nothing to scale", call. = FALSE)
    }
  }

  # deSolve needs list(dydt); a bare derivative vector is friendlier to
  # write and costs one wrapper
  func <- function(t, y, p) {
    r <- dynamics(t, y, p)
    if (is.list(r)) r else list(r)
  }

  # An infusion is a constant rate added to the derivative over one
  # segment. It rides along as `n_state` extra parameters rather than as
  # a closure constant or a forcing function: parameters are tape inputs,
  # so an estimated rate is differentiated exactly, and a constant one is
  # dropped from the sensitivity system by RTMBode itself. A branch on
  # `t` inside the dynamics would be the obvious alternative and is not
  # available - RTMBode tapes the derivative function once, so `t` is an
  # advector and RTMB refuses the comparison.
  #
  # A `tv` value rides the same idea one step earlier in the vector: the
  # dynamics is handed `c(parms, tv)`, so it reads a time-varying input
  # at `p[length(parms) + j]` and never learns that the value changed.
  n_parm <- length(parm_cols)
  n_pall <- n_parm + n_tv
  dyn_rate <- if (!has_infusion) NULL else function(t, y, p) {
    r <- dynamics(t, y, p[seq_len(n_pall)])
    if (!is.list(r)) return(list(r + p[n_pall + seq_len(n_state)]))
    r[[1L]] <- r[[1L]] + p[n_pall + seq_len(n_state)]
    r
  }

  col_at <- function(cl, i) if (length(cl) == 1L) cl else cl[i]

  # One RTMBode call plus the failure checks, shared by the single-solve
  # path and by every segment of a dosing solve.
  #
  # deSolve reports "integration was not successful" as a warning and
  # still returns a full-length matrix of whatever it reached, so the
  # warning is the only evidence that the numbers are not a solution. It
  # is left to propagate as well: during a fit it is the user's signal,
  # and there it is all we have.
  run_solve <- function(y, grid, pvals, fn) {
    gave_up <- NULL
    s <- withCallingHandlers(
      RTMBode::ode(y = y, times = grid, func = fn %||% func, parms = pvals,
                   method = method, atol = atol, rtol = rtol, ...),
      warning = function(w) {
        if (grepl(ode_giveup_pattern, conditionMessage(w))) {
          gave_up <<- conditionMessage(w)
        }
      }
    )
    # a solver that gives up can also return a short matrix; that must be
    # a failure, not a length error later
    if (nrow(s) != length(grid)) {
      stop("the integrator returned ", nrow(s), " of ", length(grid),
           " requested time points", call. = FALSE)
    }
    # The remaining two checks work on numbers only: RTMB refuses
    # comparison on AD types, so on the tape a diverging trajectory
    # cannot be seen at all. See the "Failed solves" section of the help
    # page for what that means for on_error.
    if (!inherits(s, "advector")) {
      if (!is.null(gave_up)) {
        stop("the integrator did not converge: ", gave_up, call. = FALSE)
      }
      if (!all(is.finite(s))) {
        stop("the integrator returned non-finite values", call. = FALSE)
      }
    }
    s
  }

  outs <- lapply(out_idx, function(k) numeric(n_obs))
  failed <- character(0)
  ss_acc <- new.env(parent = emptyenv())
  ss_acc$groups <- character(0)
  ss_acc$rel <- 0
  ss_acc$stalled <- FALSE

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
    y0 <- do.call(c, lapply(init_cols, col_at, i = i1))
    pv <- do.call(c, lapply(parm_cols, col_at, i = i1))
    ev <- if (is.null(ev_by_group)) NULL else ev_by_group[[labels[[g]]]]
    if (!is.null(ev) && any(ev$time < tstart)) {
      j <- which(ev$time < tstart)[1L]
      stop("group '", labels[[g]], "' has an event at time ",
           format(ev$time[j]), ", before t0 (", format(tstart),
           "); frm_ode() integrates forward from t0", call. = FALSE)
    }

    tvb <- if (!n_tv) NULL else
      ode_tv_blocks(tv_cols, brk, idx, times, tstart, labels[[g]])

    sol <- tryCatch(
      if (is.null(ev) && is.null(tvb)) {
        # deSolve requires times[1] to be the initial time; the duplicate
        # that an observation at t0 creates is tolerated, and the extra
        # row is dropped from the solution
        s <- run_solve(y0, c(tstart, times[idx]), pv, NULL)
        lapply(seq_len(n_state), function(k) s[-1L, k + 1L])
      } else {
        scale_g <- if (is.null(scale_cols)) 1 else
          col_at(scale_cols[[1L]], i1)
        ode_solve_events(run_solve, y0, pv, times[idx], ev, tstart,
                         n_state, dyn_rate, scale_g, tvb, n_ss, ss_tol,
                         ss_acc, labels[[g]], ss_extrapolate, atol,
                         rtol)
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
      z[idx] <- sol[[out_idx[k]]]
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

  # The run-in is a fixed number of cycles, so nothing forces it to have
  # arrived. Off the tape the successive extrapolants can be compared,
  # and that is the only place the shortfall is visible at all, so it is
  # said out loud there. The number reported is a distance to the LIMIT;
  # the cycle-to-cycle movement this used to print understated it by
  # about 1 / (lambda_z * ii), which is worst where the error is worst.
  if (length(ss_acc$groups)) {
    warning(
      "frm_ode(): after ", n_ss, " steady-state run-in cycles the ",
      "state at the start of the cycle is about ",
      format(signif(ss_acc$rel, 3)), " (relative) away from the limit ",
      "it is approaching, more than ss_tol = ", format(ss_tol), ", in ",
      length(ss_acc$groups), " group",
      if (length(ss_acc$groups) == 1L) "" else "s", " (",
      paste(utils::head(ss_acc$groups, 5L), collapse = ", "),
      if (length(ss_acc$groups) > 5L) ", ..." else "", "). ",
      if (isTRUE(ss_acc$stalled))
        paste0("At least one state is not contracting between cycles ",
               "at all, so it has no steady state to reach and no ",
               "n_ss is enough for it: an `ss` row assumes every state ",
               "settles. ")
      else if (ss_extrapolate)
        "Raise n_ss. "
      else
        paste0("Raise n_ss, or leave ss_extrapolate at TRUE, which ",
               "sums this tail instead of truncating it. "),
      "That figure DETECTS rather than measures: it is built from the ",
      "same geometric model the correction is, so where the model is ",
      "poor it is poor with it, and it is a lower bound rather than ",
      "an estimate. This check runs on the numeric path only - on the ",
      "tape RTMB refuses comparison, so a fit cannot repeat it",
      call. = FALSE)
  }

  if (!is.null(out_sel)) {
    # one shared solve, one state per row: start from the first state
    # everywhere and overwrite the rows that read a different one
    res <- outs[[1L]]
    for (k in seq_len(n_state)[-1L]) {
      rows <- which(out_sel == k)
      if (length(rows)) res[rows] <- outs[[k]][rows]
    }
    return(res)
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
#
# frmtmb calls this through frmtmb_register_frame_check(), which .onLoad
# fills in. The core therefore names neither this package nor frm_ode().

#' Every `frm_ode()` or `frm_lincmt()` call inside an expression,
#' matched to its formals.
#'
#' Both are covered because both read their dynamics inputs off the
#' group's first row, so both are wrong in the same silent way when a
#' covariate varies inside a group.
#'
#' @noRd
find_ode_calls <- function(expr, out = list()) {
  if (is.call(expr)) {
    fn <- expr[[1L]]
    nm <- if (is.name(fn)) as.character(fn) else
      if (is.call(fn) && identical(as.character(fn[[1L]]), "::")) {
        as.character(fn[[3L]])
      } else ""
    if (nm %in% c("frm_ode", "frm_lincmt")) {
      f <- if (identical(nm, "frm_ode")) frm_ode else frm_lincmt
      out <- c(out, list(
        tryCatch(match.call(f, expr), error = function(e) NULL)
      ))
    }
    for (i in seq_along(expr)) {
      if (i == 1L) next
      if (!missing_arg(expr[[i]])) out <- find_ode_calls(expr[[i]], out)
    }
  }
  out
}

#' Is a matched call argument the empty symbol?
#'
#' `match.call()` leaves an unsupplied argument as the empty symbol,
#' which is neither `NULL` nor missing to `is.null()`. Returns `TRUE`
#' for that value only.
#'
#' @noRd
missing_arg <- function(x) is.name(x) && !nzchar(as.character(x))

#' Refuse a dynamics parameter that varies inside a solve group.
#'
#' Walks the nonlinear body of every linear predictor for `frm_ode()`
#' and `frm_lincmt()` calls, and for each nonlinear parameter appearing
#' in that call's `init` or `parms` argument checks that the
#' parameter's design is constant within the levels of the call's
#' `group` column. Only a bare symbol `group =` can be resolved against
#' the model frame; anything computed is left to the runtime check in
#' [frm_ode()].
#'
#' The signature is `frmtmb_register_frame_check()`'s: `frame` supplies
#' the model frame and the `linpred()` accessor.
#'
#' @noRd
check_ode_constancy <- function(spec, frame) {
  mf <- frame[["data_frame"]]
  linpred <- frame[["linpred"]]
  for (resp in spec$responses) {
    # every nonlinear body of the response, not just mu's: nlf() can put
    # a solve behind any parameter
    bodies <- Filter(Negate(is.null),
                     lapply(resp$dpars, function(dp) dp[["nl_body"]]))
    calls <- unlist(lapply(bodies, function(nl) {
      Filter(Negate(is.null), find_ode_calls(nl))
    }), recursive = FALSE)
    if (!length(calls)) next
    nlpars <- resp$nlpars %||% character(0)
    for (cl in calls) {
      gexpr <- cl[["group"]]
      if (is.null(gexpr) || !is.name(gexpr)) next
      gname <- as.character(gexpr)
      gv <- mf[[gname]]
      if (is.null(gv) || is.matrix(gv)) next
      gi <- as.integer(factor(gv))
      for (arg in c("init", "parms", "t0", "event_scale")) {
        aexpr <- cl[[arg]]
        if (is.null(aexpr)) next
        for (np in intersect(all.vars(aexpr), nlpars)) {
          lp <- linpred(resp$resp_name, np)
          if (is.null(lp)) next
          bad <- c(ode_varying_cols(lp[["X"]], gi),
                   ode_varying_cols(lp[["Z"]], gi))
          if (length(bad)) {
            stop("Nonlinear parameter '", np, "' is a dynamics input of ",
                 utils::tail(as.character(cl[[1L]]), 1L),
                 "() but is not constant within '", gname,
                 "': ", paste(unique(bad), collapse = ", "),
                 ". One system is solved per group and every dynamics ",
                 "input is read off the group's first row, so a ",
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
