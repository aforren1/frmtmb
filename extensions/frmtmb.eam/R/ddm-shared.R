# Pieces more than one family in this package needs, kept in one place so
# that the second family to want them does not grow a second copy.

#' `exp(la) - exp(lb)` from the two LOGS, without overflow and without
#' cancellation.
#'
#' Both differences this package forms in a tail are differences of two
#' exponentials whose logs are known exactly: the normal density
#' difference in the racing diffusion model's density, and the
#' difference of the two boundary terms in its survival. Written
#' directly as `exp(la) - exp(lb)` either overflows, when one log is
#' large and positive, or underflows one term to zero while the other
#' stays finite, which turns a small difference into the larger term
#' alone.
#'
#' Anchoring at the larger log fixes both at once. Subtract the anchor
#' from each: one argument becomes exactly zero and the other is at most
#' zero, so neither `exp()` can overflow. Then use `expm1` on both rather
#' than `exp`, because `expm1` of the anchor's own argument is exactly
#' zero: the difference is the non-anchored term alone, computed to full
#' relative precision, rather than `1 - something-near-1`.
#'
#' `ddm_floor(la, lb)` is the anchor. It is written as a floor and used
#' here as a maximum, which is the same function: the helper is
#' `max(x, lo)` for any `lo`, and nothing in it requires `lo` to be a
#' constant.
#'
#' @noRd
ddm_expdiff <- function(la, lb) {
  m <- ddm_floor(la, lb)
  exp(m) * (expm1(la - m) - expm1(lb - m))
}

#' Log of the standard normal density, as a plain expression.
#'
#' `log(dnorm(x))` underflows to `-Inf` for `|x|` past about 38, which is
#' inside the range a race's loser reaches. The log is an exact
#' polynomial and never does.
#'
#' @noRd
ddm_lphi <- function(x) -0.5 * log(2 * pi) - 0.5 * x * x

#' The fastest response, in seconds, above which the data is read as a
#' unit mistake rather than as a slow task.
#'
#' Not a plausible single response time. It is a ceiling on the fastest
#' response IN THE WHOLE DATA SET, which is a much rarer thing: one slow
#' trial does not reach it and cannot.
#'
#' Where 20 comes from, measured on simulated wiener data
#' (`dev/smallitems-findings.md` carries the tables). Over 192 cells at
#' the parameters a two-choice task actually produces, boundary 2 to 5,
#' non-decision time 0.5 to 5 seconds, drift 0.3 to 2, at 12 to 80
#' trials, the fastest response ran from 0.599 to 6.605 seconds and
#' NONE false alarmed; and over 81 cells at the standard published
#' range, all 81 are caught when the same data is read as milliseconds.
#' In the other direction the miss rate is zero at every threshold up to
#' 100 and rises at 200, so 20 sits a factor of five inside the band
#' where nothing is missed.
#'
#' It CAN fire on a correct model, and the rate is worth knowing rather
#' than assuming. The exposed design is a slow one with a short
#' session. At 20 trials the rate is 0.045 where the fastest response
#' averages 13.4 seconds, 0.155 at 15.7, 0.690 at 22.7 and 1.000 at
#' 41.7; sixty trials of the same tasks give 0.000, 0.000, 0.285 and
#' 0.995, because the fastest of more draws settles toward the floor.
#'
#' Quoted against the FASTEST response and not the median, because the
#' median does not determine the rate: at a median near 50 seconds the
#' rate runs from 0.185 to 0.935 depending on whether the task is slow
#' from a far boundary or from weak evidence, and the fastest response
#' tracks it in every one of those cells. The cost falls on
#' deliberation, insight and matrix-reasoning designs, whose
#' non-decision time is 1 to 2 seconds rather than 20.
#'
#' @noRd
ddm_seconds_ceiling <- 20

#' Warn when the response is not on the scale the defaults assume.
#'
#' Every default in this package reads the response as SECONDS. The
#' starting values, `gddm_control(dt = 0.01)` and the window `t_max`
#' takes from the data are all absolute times, and nothing in any of the
#' five likelihoods refuses milliseconds: the fit converges and reports
#' a boundary separation three orders of magnitude out. That is the
#' silent wrong answer this package ranks first, so it is worth a
#' warning even at the cost of an occasional false alarm.
#'
#' A warning and not a refusal, because a design with no trial under 20
#' seconds is possible, if rare. It carries a class so that such a
#' design can silence this one condition without also hiding the
#' convergence warnings beside it.
#'
#' Called from each family's `valid_y` rather than once from here, so
#' that the sentence a user sees names the family they wrote.
#'
#' @noRd
ddm_check_units <- function(y, what) {
  if (!length(y)) return(invisible(NULL))
  lo <- min(y)
  if (!is.finite(lo) || lo <= ddm_seconds_ceiling) return(invisible(NULL))
  warning(warningCondition(paste0(
    what, ": the fastest response in these data is ",
    format(lo, digits = 4), ", and every default in this package reads ",
    "the response as a time in SECONDS. Milliseconds is the usual ",
    "cause: a millisecond clock puts a typical response near 500, and a ",
    "task in which no trial at all finishes within ",
    ddm_seconds_ceiling, " seconds is rare. Divide the response by 1000 ",
    "if that is what happened, which puts these times between ",
    format(lo / 1000, digits = 4), " and ",
    format(max(y) / 1000, digits = 4), " seconds. If the task really is ",
    "this slow then the fit is correct and this warning is its only ",
    "cost; it carries the class frmtmb_eam_units_warning so that it can ",
    "be silenced on its own"),
    class = "frmtmb_eam_units_warning"))
  invisible(NULL)
}

# ------------------------------------- the non-decision time and its bound
#
# WHY THERE ARE TWO PARAMETERIZATIONS. The density of every family here
# is zero at and below `ndt`, so the likelihood has a hard edge at the
# fastest response. A log link would let the optimizer walk over it; a
# logit scaled onto `(0, ub)` makes the constraint structural.
#
# One `ub` cannot express the constraint when `ndt` carries a random
# effect. The information about a subject's non-decision time is that
# subject's OWN fastest response, and at 30 subjects by 400 trials with
# a between-subject spread of 26 ms on a mean of 250 ms, 20 of the 30
# subjects have a true `ndt` above the global minimum while NONE is
# above its own; the fit then runs to the wall and reports a 7.2e-06
# standard error on an answer wrong by ten percent
# (`dev/scale-findings.md`).
#
# A per-ROW bound cannot live in a link: `linkinv()` is handed a vector
# whose length is the rows being predicted, so a captured per-row bound
# recycles silently against newdata of another length, and `linkfun()`
# is called once on a SCALAR starting value. So:
#
#   NO ndt_group():  the bound is one number and STAYS IN THE LINK, as
#     it has since 0.5.0. `ndt` is a TIME on the response scale, a
#     `prior(class = "ndt")` is a density on that time, a `bf(ndt = )`
#     constant is that time, and nothing about such a model moved.
#   WITH ndt_group(): the bound is per row, so the link becomes a plain
#     logit and `ndt` is a FRACTION of the row's own bound, which the
#     density multiplies back out. That is the breaking change, and it
#     is confined to models that opt in.
#
# The per-row bound reaches the density as DATA, `ndt_floor`, beside the
# addition-term values, so the tape gathers a vector instead of doing a
# lookup per evaluation.

#' A logit scaled onto `(0, ub)`.
#'
#' `ub` may be `NA`, which is the state a family constructed without
#' `max_ndt` is in before [frm()] hands it the data. Using such a link
#' says so rather than returning a silent `NA`, and that refusal is
#' load-bearing in two places a fit never reaches: `bf(ndt = 0.2)` runs
#' `linkfun()` at parse time, and [frmtmb::mixture()] does not finalize
#' its components at all.
#'
#' @noRd
ddm_scaled_logit <- function(ub, dpar, what = "wiener") {
  force(ub)
  force(dpar)
  force(what)
  bound <- function() {
    if (is.na(ub)) {
      stop(what, "(): the ", dpar, " bound is not set yet. This ",
           "happens when a ", what, "() family object is used outside ",
           "frm(), for example to inspect its links before a fit, to ",
           "pin the parameter with bf(", dpar, " = ), or inside ",
           "mixture(), which does not finalize its components. Pass ",
           "`max_ndt` to ", what, "() to set the bound up front.",
           call. = FALSE)
    }
    ub
  }
  list(
    # the bare word, which is what 92e9330 printed for this family and
    # what two tests pin. The bound itself is on the fitted family, at
    # `family(fit)$ndt_bound$ub`, which is where ?wiener sends a reader
    # rather than into a link's display string.
    name = "scaled_logit",
    linkfun = function(mu) { U <- bound(); log(mu / (U - mu)) },
    linkinv = function(eta) bound() / (1 + exp(-eta)),
    mu_eta = function(eta) {
      p <- 1 / (1 + exp(-eta))
      bound() * p * (1 - p)
    })
}

#' A code for a group LABEL that does not depend on which other labels
#' are present.
#'
#' The first version of this coerced a factor with `as.integer()`, that
#' is by the level INDEX, arguing that a data frame carries its levels
#' through a subset. It carries the level SET through `[`. It does not
#' carry the ORDER through `droplevels()`, `relevel()`, or a prediction
#' grid the user builds, and each of those permutes the indices. Review
#' `dev/reviews/2026-09-09-ndt.md` B1 measured it: after `droplevels()`
#' every row was paired with the PREVIOUS group's bound, silently, which
#' is the defect this whole change exists to remove.
#'
#' So the code is a function of the label and of nothing else: two
#' polynomial rolling hashes over the label's code points, packed into
#' one double below 2^49 so the arithmetic is exact. Collisions inside
#' one column are refused where they can be seen, in `ddm_ndt_group()`,
#' which has the labels. A collision BETWEEN a fitted label and a label
#' seen only on newdata cannot be seen from here and would pair that row
#' with the collided group instead of refusing it as unseen; at 2^-49
#' per unseen label that is far below every other error in this package,
#' and it is stated rather than hidden.
#'
#' @noRd
ddm_label_code <- function(lab) {
  vapply(lab, function(s) {
    cp <- utf8ToInt(s)
    if (!length(cp)) cp <- 0L
    h1 <- 0
    h2 <- 1
    for (k in cp) {
      h1 <- (h1 * 131 + k) %% 33554393
      h2 <- (h2 * 137 + k) %% 16777213
    }
    h1 * 16777216 + h2
  }, numeric(1), USE.NAMES = FALSE)
}

#' Coerce an `ndt_group()` column to those codes.
#'
#' Every spelling goes through the LABEL, so a factor, a character
#' vector, a logical and an integer column that name the same groups all
#' give the same fit, and none of them depends on the level order of the
#' frame it was evaluated in.
#'
#' @noRd
ddm_coerce_ndt_group <- function(x) {
  lab <- if (is.factor(x)) as.character(x) else as.character(x)
  if (anyNA(lab)) {
    stop("ndt_group(): the grouping says which trials share a ",
         "non-decision-time bound, so every row needs one, and this ",
         "column has a missing value. frm()'s own na.action drops such ",
         "a row before this is reached; predicting on newdata does ",
         "not, which is where this fires.", call. = FALSE)
  }
  u <- unique(lab)
  co <- ddm_label_code(u)
  if (anyDuplicated(co)) {
    j <- which(duplicated(co) | duplicated(co, fromLast = TRUE))
    stop("ndt_group(): the group labels ",
         paste(sort(u[j])[seq_len(min(4L, length(j)))], collapse = ", "),
         " collide under the code this term is keyed on, so they cannot ",
         "be told apart. Rename one of them.", call. = FALSE)
  }
  co[match(lab, u)]
}

#' The bound each row's non-decision time is measured against.
#'
#' Returns the scalar bound and, when `ndt_group()` is present, the
#' per-group table of fastest responses keyed by the label codes.
#' Refusing the pair is the whole of `max_ndt`'s meaning under the new
#' scheme: `max_ndt` sets the bound to one number the user chose and
#' `ndt_group()` sets it to a number per group that the data chose, and
#' a model cannot want both.
#'
#' @noRd
ddm_ndt_spec <- function(y, aterms, max_ndt, what) {
  g <- aterms[["ndt_group"]]
  if (!is.null(g) && !is.null(max_ndt)) {
    stop(what, ": max_ndt and ndt_group() both set the non-decision ",
         "time's upper bound and they set it to different things. ",
         "max_ndt fixes one bound for every row; ndt_group() takes ",
         "each group's own fastest response. Drop one of them.",
         call. = FALSE)
  }
  if (is.null(g)) return(list(ub = max_ndt %||% min(y), floors = NULL))
  # One value is allowed and means one group, which is the same model as
  # no ndt_group() at all. Anything else that is not one per trial is a
  # recycling accident rather than a design.
  if (!length(g) %in% c(1L, length(y))) {
    stop(what, ": ndt_group() has ", length(g), " values for ",
         length(y), " responses. It names the group of every trial.",
         call. = FALSE)
  }
  k <- as.character(g)
  list(ub = min(y),
       floors = vapply(split(as.numeric(y), k), min, numeric(1)),
       sizes = lengths(split(as.numeric(y), k)))
}

#' The per-row bound, wherever the row came from.
#'
#' Only a GROUPED model has one of these. `ndt_floor` is the vector
#' frame assembly baked in and is what the tape reads; failing that the
#' row's group is looked up in the table this closure captured at
#' finalize, which is the newdata path. The bound is a property of the
#' FITTED data, so a group the fit never saw has no bound and is refused
#' rather than given the global one.
#'
#' @noRd
ddm_ndt_scaler <- function(floors, what) {
  force(floors)
  force(what)
  function(aterms) {
    fl <- aterms[["ndt_floor"]]
    if (!is.null(fl)) return(fl)
    g <- aterms[["ndt_group"]]
    if (is.null(g)) {
      stop(what, ": this model bounds the non-decision time by each ",
           "ndt_group()'s own fastest response, so a prediction needs ",
           "every row's group. Supply the ndt_group() column on ",
           "newdata.", call. = FALSE)
    }
    v <- floors[as.character(g)]
    if (anyNA(v)) {
      stop(what, ": ", sum(is.na(v)), " row(s) name an ndt_group() the ",
           "model was not fitted to, so there is no fastest response to ",
           "bound their non-decision time by. A per-group bound is a ",
           "property of the fitted data. For wiener_gng() a group whose ",
           "trials are all no-go rows lands here too, because the bound ",
           "is taken over the go rows alone.", call. = FALSE)
    }
    unname(v)
  }
}

#' Put `ndt` (and `st`) back on the response's own scale.
#'
#' GROUPED models only. Every slot that consumes a non-decision time
#' goes through here first, so each density, mean and simulator below is
#' written in TIMES and knows nothing about the fraction. `st` is a
#' width on the same quantity and is a fraction of twice the bound,
#' which is the range its own scaled logit carries in the scalar case.
#'
#' The reserved `.eta_` entries are dropped with the value they belong
#' to. They are the linear predictor the ROBUST dpar accessors recover a
#' saturated value from, and after this the linear predictor beside
#' `ndt` no longer maps to it; leaving a stale one there would let an
#' accessor return a different number from the one the density used.
#'
#' @noRd
ddm_ndt_natural <- function(dpars, aterms, scale) {
  fl <- scale(aterms)
  dpars[["ndt"]] <- dpars[["ndt"]] * fl
  dpars[[".eta_ndt"]] <- NULL
  if (!is.null(dpars[["st"]])) {
    dpars[["st"]] <- dpars[["st"]] * (2 * fl)
    dpars[[".eta_st"]] <- NULL
  }
  dpars
}

#' Wrap one slot so that it receives times rather than fractions.
#'
#' The formal COUNT is preserved rather than replaced by `...`, because
#' frmtmb reads `length(formals())` to decide whether a slot takes the
#' family-level extras (`R/families.R`, `R/objective.R`). A wrapper with
#' a dots formal would have been read as taking one argument.
#'
#' @noRd
ddm_wrap_q <- function(f, scale) {
  force(f)
  force(scale)
  if (length(formals(f)) >= 4L) {
    function(y, dpars, aterms, extra) {
      f(y, ddm_ndt_natural(dpars, aterms, scale), aterms, extra)
    }
  } else {
    function(y, dpars, aterms) {
      f(y, ddm_ndt_natural(dpars, aterms, scale), aterms)
    }
  }
}

#' @noRd
ddm_wrap_sim <- function(f, scale) {
  force(f)
  force(scale)
  if (length(formals(f)) >= 4L) {
    function(dpars, aterms, n, extra) {
      f(ddm_ndt_natural(dpars, aterms, scale), aterms, n, extra)
    }
  } else {
    function(dpars, aterms, n) {
      f(ddm_ndt_natural(dpars, aterms, scale), aterms, n)
    }
  }
}

#' @noRd
ddm_wrap_mean <- function(f, scale) {
  force(f)
  force(scale)
  function(dpars, aterms) {
    f(ddm_ndt_natural(dpars, aterms, scale), aterms)
  }
}

#' Install the bound on a family that is otherwise finished.
#'
#' Three states, and only the third is new arithmetic.
#'
#' - PENDING, no bound yet: the scaled logit refuses, which is what a
#'   family object used outside `frm()` did before this change and what
#'   `bf(ndt = )` and `mixture()` still meet.
#' - SCALAR bound: the scaled logit carries it, so the density, the
#'   starting values, a prior on `class = "ndt"` and a pinned constant
#'   all see a TIME. Byte for byte the pre-change path.
#' - PER-GROUP bound: a plain logit on the fraction, the slots wrapped,
#'   and `ndt_floor` added to the addition-term values.
#'
#' Idempotent, and that is not decoration: `assemble_frame()` runs
#' `family_finalize()` again on the family it returned last time, on the
#' `influence()`, `frm_simulate()` and prior-predictive paths, and a
#' second wrap would multiply the non-decision time by its bound twice.
#' A settled bound is also KEPT rather than re-derived, so a refit on
#' fewer rows scores against the bound the fit was made with; the four
#' families reach this through `ddm_ndt_keep()`.
#'
#' @noRd
ddm_ndt_install <- function(fam, ub, floors, what, pending = FALSE,
                            y_of = ddm_all_rows, sizes = NULL) {
  raw <- fam[["ndt_raw"]]
  if (is.null(raw)) {
    raw <- list(lpdf = fam[["lpdf"]], lcdf = fam[["lcdf"]],
                lccdf = fam[["lccdf"]], sim = fam[["sim"]],
                mean_fn = fam[["post"]][["mean_fn"]],
                aterm_data = fam[["aterm_data"]],
                init_ndt = fam[["init_dpars"]][["ndt"]],
                init_st = fam[["init_dpars"]][["st"]])
  } else if (!isTRUE(fam[["ndt_bound"]][["pending"]])) {
    return(fam)   # already carrying a settled bound
  }
  has_st <- !is.null(raw[["init_st"]])
  # start every slot from the UNWRAPPED originals, so a re-install over
  # a pending bound cannot stack a second wrapper
  for (nm in c("lpdf", "lcdf", "lccdf")) fam[[nm]] <- raw[[nm]]
  fam[["sim"]] <- raw[["sim"]]
  fam[["post"]][["mean_fn"]] <- raw[["mean_fn"]]
  fam[["aterm_data"]] <- raw[["aterm_data"]]

  if (is.null(floors)) {
    # the scalar case, which is every model that predates ndt_group():
    # the bound goes back in the link and nothing else moves
    lo <- if (pending) NA_real_ else ub
    fam[["links"]][["ndt"]] <- ddm_scaled_logit(lo, "non-decision time",
                                                what)
    if (has_st) {
      fam[["links"]][["st"]] <-
        ddm_scaled_logit(if (is.na(lo)) NA_real_ else 2 * lo,
                         "non-decision time range", what)
    }
    fam[["init_dpars"]][["ndt"]] <- ddm_ndt_init_time(y_of)
    if (has_st) fam[["init_dpars"]][["st"]] <- ddm_st_init_time(y_of)
  } else {
    scale <- ddm_ndt_scaler(floors, what)
    fam[["links"]][["ndt"]] <- "logit"
    if (has_st) fam[["links"]][["st"]] <- "logit"
    for (nm in c("lpdf", "lcdf", "lccdf")) {
      if (is.function(raw[[nm]])) fam[[nm]] <- ddm_wrap_q(raw[[nm]], scale)
    }
    if (is.function(raw[["sim"]])) {
      fam[["sim"]] <- ddm_wrap_sim(raw[["sim"]], scale)
    }
    if (is.function(raw[["mean_fn"]])) {
      fam[["post"]][["mean_fn"]] <- ddm_wrap_mean(raw[["mean_fn"]], scale)
    }
    # The floor is data and rides with the addition terms, so the tape
    # gathers it once instead of looking a group up per evaluation. It
    # is read out of the CAPTURED table rather than recomputed from the
    # `y` this is called with.
    prev <- raw[["aterm_data"]]
    fam[["aterm_data"]] <- if (is.null(prev)) {
      function(y, aterms) list(ndt_floor = scale(aterms))
    } else {
      function(y, aterms) {
        c(prev(y, aterms), list(ndt_floor = scale(aterms)))
      }
    }
    # half of each group's own bound, which is what 0.5 * min(y) meant
    # when there was one bound
    fam[["init_dpars"]][["ndt"]] <- function(y, aterms) 0.5
    if (has_st) fam[["init_dpars"]][["st"]] <- function(y, aterms) 0.05
  }
  fam[["ndt_raw"]] <- raw
  fam[["ndt_bound"]] <- list(ub = ub, floors = floors, pending = pending,
                             sizes = sizes)
  fam
}

#' Half the fastest response, in the response's own units.
#'
#' The starting value the scaled logit has always taken. Expressed
#' against `min(y)` rather than against the bound, because under a
#' `max_ndt` ABOVE the fastest response, which is the mixture case, half
#' the BOUND starts the optimizer at a non-decision time no observed row
#' can reach: measured on `test-defects.R`'s mixture, that start walked
#' the fit to `ndt` = 2.7e-07 with a drift of 51.7, a HIGHER log
#' likelihood at a degenerate point.
#'
#' @noRd
ddm_ndt_init_time <- function(y_of) {
  force(y_of)
  function(y, aterms) 0.5 * min(y_of(y, aterms))
}

#' @noRd
ddm_st_init_time <- function(y_of) {
  force(y_of)
  function(y, aterms) 0.1 * min(y_of(y, aterms))
}

#' Which rows a starting value and a bound may read. Every family but
#' the go/no-go one reads all of them; that one reads the GO rows,
#' because a no-go row's response entry is a placeholder and letting one
#' into a minimum would let it set a parameter's range.
#'
#' @noRd
ddm_all_rows <- function(y, aterms) y

#' The bound a family carries before [frm()] has seen a response.
#'
#' `max_ndt` settles it up front and nothing else can, so a family built
#' without one carries a refusing link. That refusal is not decoration:
#' `bf(ndt = 0.2)` runs `linkfun()` at parse time and
#' [frmtmb::mixture()] never finalizes its components, so both would
#' otherwise read a number on a scale nothing had set.
#'
#' @noRd
ddm_ndt_preinstall <- function(fam, max_ndt, what,
                               y_of = ddm_all_rows) {
  if (is.null(max_ndt)) {
    return(ddm_ndt_install(fam, NA_real_, NULL, what, pending = TRUE,
                           y_of = y_of))
  }
  ddm_ndt_install(fam, max_ndt, NULL, what, y_of = y_of)
}

#' The bound a family is already carrying, or `NULL`.
#'
#' A family whose bound is settled KEEPS it: `assemble_frame()` re-runs
#' `family_finalize()` on the `influence()`, `frm_simulate()` and
#' prior-predictive paths, and a bound re-derived from fewer rows would
#' make a leave-one-out refit a refit of a different model. Families
#' that rebuild themselves from their config have to ask for this
#' explicitly, which is why it is a function rather than a side effect
#' of the install.
#'
#' @noRd
ddm_ndt_keep <- function(fam) {
  bd <- fam[["ndt_bound"]]
  if (is.null(bd) || isTRUE(bd[["pending"]])) return(NULL)
  bd
}

#' Fit a bounded non-decision-time link to the observed response.
#'
#' Shared rather than copied because four families in this package want
#' the same bound and the same refusal, and the refusal carries the
#' family's name at run time so that the message a user sees still names
#' the family they wrote.
#'
#' @noRd
ddm_ndt_finalize <- function(fam, y, aterms, max_ndt, what,
                             y_of = ddm_all_rows) {
  sp <- ddm_ndt_spec(y, aterms, max_ndt, what)
  if (!is.null(max_ndt) && sp$ub > min(y)) {
    stop(what, ": max_ndt = ", format(sp$ub), " is above the fastest ",
         "response (", format(min(y)), "). Nothing can be observed ",
         "before the non-decision time, so a bound above the fastest ",
         "response admits values at which that trial has no ",
         "likelihood.", call. = FALSE)
  }
  keep <- ddm_ndt_keep(fam)
  if (!is.null(keep)) return(fam)
  ddm_ndt_install(fam, sp$ub, sp$floors, what, y_of = y_of,
                  sizes = sp$sizes)
}

#' The non-decision time, in the response's own units
#'
#' Without `ndt_group()` this is `predict(fit, dpar = "ndt", type =
#' "response")`, which already reports a time. With `ndt_group()` the
#' `ndt` link is a plain logit on a FRACTION of the row's own bound, so
#' `predict()` reports that fraction and this multiplies it back out.
#' Written so that a caller need not know which of the two a fit is.
#'
#' The bound is the fastest response of the row's `ndt_group()`, or of
#' the whole data set when the model has no `ndt_group()`, or the
#' `max_ndt` the family was given. Whichever it is, it is a property of
#' the data the model was FITTED to, so a prediction on new rows is
#' scaled by the same bound the fit used.
#'
#' @param object A fitted [wiener()], [lba()], [rdm()] or [wiener_gng()]
#'   model.
#' @param newdata Optional data frame. It must carry the model's
#'   `ndt_group()` column when the model has one.
#' @param ... Passed to [stats::predict()]; `re.form` and
#'   `allow_new_levels` are the useful ones.
#'
#' @return A numeric vector, one non-decision time per row.
#' @seealso [wiener()] for what the two parameterizations are and why
#'   there are two.
#' @examples
#' set.seed(1)
#' d <- ddm_simulate(300, mu = 1.2, bs = 1.5, ndt = 0.25)
#' fit <- frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
#'            family = wiener(), data = d)
#' # with no ndt_group() the two agree: `ndt` is already a time
#' head(predict(fit, dpar = "ndt", type = "response"), 3)
#' head(ndt_time(fit), 3)
#' @export
ndt_time <- function(object, newdata = NULL, ...) {
  rsp <- frmtmb::single_response(object)
  fam <- rsp[["family"]]
  bd <- fam[["ndt_bound"]]
  if (is.null(bd)) {
    stop("ndt_time(): a ", fam[["family"]], " model does not estimate ",
         "the non-decision time against a bound this can rescale. ",
         "ndt_time() reports for wiener(), lba(), rdm() and ",
         "wiener_gng(); on a gddm() fit predict(dpar = \"ndt\", type = ",
         "\"response\") already gives the time.", call. = FALSE)
  }
  out <- stats::predict(object, newdata = newdata, dpar = "ndt",
                        type = "response", ...)
  if (is.null(bd[["floors"]])) return(as.numeric(out))
  scale <- ddm_ndt_scaler(bd[["floors"]], "ndt_time()")
  fl <- if (is.null(newdata)) {
    scale(object$frame[["aterm_values"]][[rsp[["resp_name"]]]])
  } else {
    ex <- rsp[["aterms"]][["ndt_group"]]
    g <- tryCatch(eval(ex, newdata, rsp[["formula_env"]]),
                  error = function(e) NULL)
    if (is.null(g)) {
      stop("ndt_time(): this model bounds the non-decision time by each ",
           "ndt_group()'s own fastest response, and ",
           deparse1(ex), " could not be evaluated on newdata. Supply ",
           "that column.", call. = FALSE)
    }
    scale(list(ndt_group = ddm_coerce_ndt_group(g)))
  }
  # An in-sample prediction comes back padded for the rows `na.action`
  # dropped and the frame's own vectors do not, so the two are aligned
  # on the padding rather than recycled against it.
  if (length(fl) > 1L && length(fl) != length(out)) {
    keep <- !is.na(out)
    if (sum(keep) != length(fl)) {
      stop("ndt_time(): ", length(fl), " bounds for ", sum(keep),
           " predicted rows. Supply newdata explicitly.", call. = FALSE)
    }
    pad <- rep(NA_real_, length(out))
    pad[keep] <- fl
    fl <- pad
  }
  as.numeric(out) * as.numeric(fl)
}

#' What each family in this package reads, in one place.
#'
#' The constructor passes its entry to `accepts_aterms` and the
#' compatibility rows derive the refusals from the same entry, so the
#' table cannot promise a term frame assembly refuses, nor refuse one it
#' takes. Written here rather than read back off the family objects
#' because the rules builder runs inside `.onLoad()`, where building
#' five families to ask them one question each is work nobody needs on a
#' package load.
#'
#' Spelled in TERMS, without parentheses, which is the vocabulary
#' `frmtmb_family(accepts_aterms =)` uses.
#'
#' @noRd
ddm_accepts <- list(
  wiener     = c("dec", "vint", "weights", "ndt_group"),
  # gddm has no ndt_group entry, and the exclusion is on SCOPE. Its
  # solver reads every dpar at the FIRST ROW of a condition
  # (`gd_densities()`), and `?gddm` already requires every row sharing a
  # condition to share every parameter value, so a valid model whose
  # `ndt` varies by subject already has a condition per subject and a
  # per-CONDITION bound would reach the density exactly as `ndt` does.
  # What that costs is a Fokker-Planck solve per subject. Not that a
  # bound could not reach it. See dev/ndt-findings.md.
  gddm       = c("dec", "vint", "vreal", "weights"),
  lba        = c("vint", "weights", "ndt_group"),
  rdm        = c("vint", "weights", "cens", "trunc", "ndt_group"),
  wiener_gng = c("dec", "vreal", "weights", "cens", "ndt_group"))
