# brms-style prior specification for frm_sample() (and, later, MAP).

#' Set up priors brms-style
#'
#' Builds prior specifications with brms spelling:
#' `set_prior("normal(0, 5)", class = "b")`. Combine several with `+` or
#' `c()`. Distributions: `normal(mu, sd)`, `student_t(df, mu, sd)`,
#' `cauchy(mu, sd)`, `exponential(rate)`; an empty string sets bounds
#' only.
#'
#' Classes and their scales:
#' - `"b"`: population-level coefficients of `dpar` (default: the
#'   location parameters), excluding the intercept; narrow to one
#'   coefficient with `coef`. Link scale.
#' - `"Intercept"`: the intercept of `dpar`. Link scale.
#' - `"sd"`: random-effect standard deviations (and smoothing SDs), on
#'   the NATURAL sd scale with the log-Jacobian applied, so
#'   `set_prior("exponential(1)", class = "sd")` means what it says;
#'   narrow with `group`. Correlation parameters are not covered
#'   (use class `"theta"` on the internal scale if you must).
#' - `"theta"`: raw internal covariance parameters (escape hatch).
#'
#' When priors overlap, later specifications override earlier ones, so
#' put class-wide priors first and coefficient-specific ones after.
#' `lb`/`ub` become hard bounds (via Stan's constrained transforms in
#' [frm_sample()]); for class `"sd"` they apply on the sd scale.
#'
#' @param prior Distribution string, e.g. `"normal(0, 5)"`, or a
#'   [prior_normal()]/[prior_t()] object, or `""` for bounds only.
#' @param class `"b"`, `"Intercept"`, `"sd"`, or `"theta"`.
#' @param coef Restrict to one coefficient (classes `"b"`/`"Intercept"`).
#' @param dpar Distributional parameter (default: the location
#'   parameters).
#' @param group Restrict class `"sd"` to one grouping factor.
#' @param lb,ub Optional hard bounds.
#' @return A `frmtmb_priorlist`.
#'
#' @srrstats {G2.0,G2.1} `prior` is asserted to be a length-one character
#'   vector before it is parsed, and the parsed distribution's arguments
#'   are asserted to have the arity that distribution requires (two for
#'   `normal`, three for `student_t`, one for `exponential`). A call that
#'   supplies neither a distribution nor bounds errors instead of
#'   producing an empty prior.
#' @srrstats {G2.3a} `class` is restricted with `match.arg()` to
#'   `"b"`, `"Intercept"`, `"sd"`, and `"theta"`, so an unexpected class
#'   errors and names the permitted values.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), z = rnorm(100),
#'                  g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'
#' # `+` combines specifications; the class-wide one goes first so the
#' # coefficient-specific one can override it
#' pr <- set_prior("normal(0, 1)", class = "b") +
#'   set_prior("normal(0, 0.2)", class = "b", coef = "z") +
#'   set_prior("exponential(1)", class = "sd", group = "g")
#' pr
#'
#' # the priors penalize the likelihood: the fit is a MAP estimate
#' fit <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, priors = pr)
#' fixef(fit)$mu
#' # the tight prior on z shrinks it toward zero
#' fixef(frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd))$mu
#'
#' # an empty distribution string sets a hard bound only
#' set_prior("", class = "b", coef = "x", lb = 0)
#'
#' # get_prior() shows which rows a design offers
#' get_prior(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
#' @export
set_prior <- function(prior = "", class = "b", coef = "", dpar = "",
                      group = "", lb = NA, ub = NA) {
  dist <- parse_prior_dist(prior)
  if (is.null(dist) && is.na(lb) && is.na(ub)) {
    stop("set_prior() needs a distribution, bounds, or both",
         call. = FALSE)
  }
  class <- match.arg(class, c("b", "Intercept", "sd", "theta"))
  spec <- list(dist = dist, class = class, coef = coef, dpar = dpar,
               group = group, lb = lb, ub = ub)
  structure(list(spec), class = "frmtmb_priorlist")
}

parse_prior_dist <- function(prior) {
  if (inherits(prior, "frmtmb_prior")) return(prior)
  stopifnot(is.character(prior), length(prior) == 1)
  if (prior == "") return(NULL)
  m <- regmatches(prior,
                  regexec("^\\s*([a-z_]+)\\s*\\(([^)]*)\\)\\s*$", prior))[[1]]
  if (length(m) != 3) {
    stop("Cannot parse prior '", prior,
         "'; expected e.g. \"normal(0, 5)\"", call. = FALSE)
  }
  kind <- m[2]
  pars <- as.numeric(strsplit(m[3], ",")[[1]])
  if (anyNA(pars)) {
    stop("Non-numeric arguments in prior '", prior, "'", call. = FALSE)
  }
  switch(kind,
    normal = {
      stopifnot(length(pars) == 2)
      prior_normal(pars[1], pars[2])
    },
    student_t = {
      stopifnot(length(pars) == 3)
      prior_t(pars[1], pars[2], pars[3])
    },
    cauchy = {
      stopifnot(length(pars) == 2)
      prior_t(1, pars[1], pars[2])
    },
    exponential = {
      stopifnot(length(pars) == 1, pars[1] > 0)
      structure(list(kind = "exponential", rate = pars[1]),
                class = "frmtmb_prior")
    },
    stop("Unsupported prior distribution '", kind,
         "' (supported: normal, student_t, cauchy, exponential)",
         call. = FALSE)
  )
}

#' @export
"+.frmtmb_priorlist" <- function(e1, e2) {
  stopifnot(inherits(e1, "frmtmb_priorlist"),
            inherits(e2, "frmtmb_priorlist"))
  structure(c(unclass(e1), unclass(e2)), class = "frmtmb_priorlist")
}

#' @export
c.frmtmb_priorlist <- function(...) {
  structure(do.call(c, lapply(list(...), unclass)),
            class = "frmtmb_priorlist")
}

#' @export
print.frmtmb_priorlist <- function(x, ...) {
  for (s in unclass(x)) {
    d <- if (is.null(s$dist)) "(bounds only)" else
      paste0(s$dist$kind, "(",
             paste(unlist(s$dist[-1]), collapse = ", "), ")")
    cat(d, " class=", s$class,
        if (nzchar(s$coef)) paste0(" coef=", s$coef),
        if (nzchar(s$dpar)) paste0(" dpar=", s$dpar),
        if (nzchar(s$group)) paste0(" group=", s$group),
        if (!is.na(s$lb)) paste0(" lb=", s$lb),
        if (!is.na(s$ub)) paste0(" ub=", s$ub), "\n", sep = "")
  }
  invisible(x)
}

#' Enumerate the targetable prior slots
#'
#' The [set_prior()] counterpart of brms's `get_prior()`: one row per
#' slot a prior can target, with the class/coef/dpar/group values to
#' pass to `set_prior()`. The default in every slot is flat (this is
#' maximum likelihood until priors are set). Class `"sd"` is targeted
#' by `group` only; class `"theta"` rows name the raw internal
#' covariance parameters (escape hatch, including correlations).
#'
#' @param formula A `bf()` formula (with family), a plain formula, or
#'   an already fitted `frmtmb_fit`.
#' @param data Model data (ignored when `formula` is a fit).
#' @param family Family, when `formula` does not carry one.
#' @return A data frame with columns `prior`, `class`, `coef`,
#'   `group`, `dpar`, `resp`, `lb`, `ub`.
#' @examples
#' dd <- data.frame(y = rnorm(60), x = rnorm(60),
#'                  g = factor(rep(1:6, 10)))
#' get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' @export
get_prior <- function(formula, data = NULL, family = NULL) {
  if (inherits(formula, "frmtmb_fit")) {
    spec <- formula$spec
    frame <- formula$frame
  } else {
    bform <- as_bform(formula, family)
    spec <- parse_spec(bform)
    frame <- assemble_frame(spec, data)
  }

  multi <- length(spec$responses) > 1L
  rows <- list()
  add <- function(class, coef = "", group = "", dpar = "", resp = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      prior = "(flat)", class = class, coef = coef, group = group,
      dpar = dpar, resp = resp, lb = NA_real_, ub = NA_real_
    )
  }

  for (lp in frame$linpreds) {
    if (!is.null(lp$constant) || !is.null(lp$nl_body)) next
    rspec <- spec$responses[[lp$resp]]
    # location dpars are the default target (dpar = ""), matching
    # set_prior()'s resolution
    dpar_lab <- if (lp$dpar %in% rspec$primary_dpars) "" else lp$dpar
    resp_lab <- if (multi) lp$resp else ""
    cn <- colnames(lp$X)
    if ("(Intercept)" %in% cn) {
      add("Intercept", dpar = dpar_lab, resp = resp_lab)
    }
    others <- setdiff(cn, "(Intercept)")
    if (length(others)) {
      add("b", dpar = dpar_lab, resp = resp_lab)
      for (co in others) {
        add("b", coef = co, dpar = dpar_lab, resp = resp_lab)
      }
    }
  }

  sd_groups <- character(0)
  for (bk in frame$re_blocks) {
    sd_i <- covstruct_registry[[bk$covstruct]]$sd_idx(bk$dim)
    if (length(sd_i)) sd_groups <- c(sd_groups, bk$group_name)
  }
  if (length(sd_groups)) {
    add("sd")
    for (g in unique(sd_groups)) add("sd", group = g)
  }
  n_th <- length(frame$par_template$theta %||% numeric(0))
  if (n_th) {
    add("theta")
    for (i in seq_len(n_th)) add("theta", coef = paste0("theta_", i))
  }

  out <- unique(do.call(rbind, rows))
  rownames(out) <- NULL
  out
}

# Resolve a priorlist against a fit: per-parameter prior entries (later
# specs override earlier) plus named bound vectors.
# Entry: comp, idx (scalar), dist, scale ("internal" or "sd"), and
# lb/ub on the ENTRY's scale (frm_simulate() rejects draws outside
# them; frm_sample() uses the named bound vectors instead).
# A later bounds-only specification tightens an entry the distribution
# already created, instead of being lost.
entry_bounds <- function(entry, s) {
  if (!is.na(s$lb)) entry$lb <- s$lb
  if (!is.na(s$ub)) entry$ub <- s$ub
  entry
}

resolve_priorlist <- function(fit, pl) {
  frame <- fit$frame
  assigned <- list()   # key "comp.idx" -> entry
  lower <- c()
  upper <- c()
  nm_of <- function(comp, idx) paste0(comp, ".", idx)

  target_coefs <- function(s) {
    # (comp, idx, name) triplets for classes b / Intercept
    out <- list()
    for (lp in frame$linpreds) {
      if (!is.null(lp$constant) || !is.null(lp$nl_body)) next
      rspec <- fit$spec$responses[[lp$resp]]
      is_loc <- lp$dpar %in% rspec$primary_dpars
      if (nzchar(s$dpar)) {
        if (!identical(lp$dpar, s$dpar)) next
      } else if (!is_loc) next
      cn <- colnames(lp$X)
      pick <- if (s$class == "Intercept") {
        which(cn == "(Intercept)")
      } else if (nzchar(s$coef)) {
        which(cn == s$coef)
      } else {
        which(cn != "(Intercept)")
      }
      for (k in pick) {
        out[[length(out) + 1L]] <- list(comp = lp$par,
                                        idx = lp$idx[k], name = cn[k])
      }
    }
    if (!length(out) && (nzchar(s$coef) || s$class == "Intercept")) {
      stop("Prior target not found (class=", s$class,
           if (nzchar(s$coef)) paste0(", coef=", s$coef),
           if (nzchar(s$dpar)) paste0(", dpar=", s$dpar), ")",
           call. = FALSE)
    }
    out
  }

  for (s in unclass(pl)) {
    if (s$class %in% c("b", "Intercept")) {
      for (tg in target_coefs(s)) {
        key <- nm_of(tg$comp, tg$idx)
        if (!is.null(s$dist)) {
          assigned[[key]] <- list(comp = tg$comp, idx = tg$idx,
                                  dist = s$dist, scale = "internal",
                                  lb = s$lb, ub = s$ub)
        } else if (!is.null(assigned[[key]])) {
          assigned[[key]] <- entry_bounds(assigned[[key]], s)
        }
        if (!is.na(s$lb)) lower[tg$name] <- s$lb
        if (!is.na(s$ub)) upper[tg$name] <- s$ub
      }
    } else if (s$class == "sd") {
      hit <- FALSE
      for (bk in frame$re_blocks) {
        if (nzchar(s$group) && !identical(bk$group_name, s$group)) next
        sd_i <- covstruct_registry[[bk$covstruct]]$sd_idx(bk$dim)
        for (k in sd_i) {
          hit <- TRUE
          i <- bk$theta_idx[k]
          key <- nm_of("theta", i)
          if (!is.null(s$dist)) {
            assigned[[key]] <- list(comp = "theta", idx = i,
                                    dist = s$dist, scale = "sd",
                                    lb = s$lb, ub = s$ub)
          } else if (!is.null(assigned[[key]])) {
            assigned[[key]] <- entry_bounds(assigned[[key]], s)
          }
          nm_theta <- paste0("theta_", i)
          if (!is.na(s$lb)) {
            lower[nm_theta] <- if (s$lb > 0) log(s$lb) else -Inf
          }
          if (!is.na(s$ub)) upper[nm_theta] <- log(s$ub)
        }
      }
      if (!hit) {
        stop("No random-effect SDs match class='sd'",
             if (nzchar(s$group)) paste0(", group=", s$group),
             call. = FALSE)
      }
    } else if (s$class == "theta") {
      n_th <- length(frame$par_template$theta %||% numeric(0))
      idx <- if (nzchar(s$coef)) {
        as.integer(sub("^theta_", "", s$coef))
      } else {
        seq_len(n_th)
      }
      for (i in idx) {
        key <- nm_of("theta", i)
        if (!is.null(s$dist)) {
          assigned[[key]] <-
            list(comp = "theta", idx = i, dist = s$dist,
                 scale = "internal", lb = s$lb, ub = s$ub)
        } else if (!is.null(assigned[[key]])) {
          assigned[[key]] <- entry_bounds(assigned[[key]], s)
        }
        if (!is.na(s$lb)) lower[paste0("theta_", i)] <- s$lb
        if (!is.na(s$ub)) upper[paste0("theta_", i)] <- s$ub
      }
    }
  }
  list(entries = unname(assigned), lower = lower, upper = upper)
}

# Log density of one prior entry value (AD-safe), with the sd-scale
# change of variables where requested.
prior_logdens <- function(x, dist, scale) {
  jac <- 0
  if (identical(scale, "sd")) {
    jac <- x          # theta = log sd; add the Jacobian
    x <- exp(x)
  }
  ld <- switch(dist$kind,
    normal = RTMB::dnorm(x, dist$location, dist$scale, log = TRUE),
    t = RTMB::dt((x - dist$location) / dist$scale, df = dist$df,
                 log = TRUE) - log(dist$scale),
    exponential = log(dist$rate) - dist$rate * x
  )
  ld + jac
}

# Accepts the legacy named list of prior objects OR a priorlist; returns
# entries + bounds.
resolve_prior_input <- function(fit, priors) {
  if (inherits(priors, "frmtmb_priorlist")) {
    return(resolve_priorlist(fit, priors))
  }
  legacy <- resolve_priors(fit, priors)
  entries <- list()
  for (e in legacy) {
    for (i in e$idx) {
      entries[[length(entries) + 1L]] <-
        list(comp = e$comp, idx = i, dist = e$prior, scale = "internal")
    }
  }
  list(entries = entries, lower = c(), upper = c())
}