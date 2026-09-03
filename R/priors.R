# brms-style prior specification for frm_sample() (and, later, MAP).

#' Set up priors brms-style
#'
#' Builds prior specifications with brms spelling:
#' `set_prior("normal(0, 5)", class = "b")`. Combine several with `+` or
#' `c()`. Distributions: `normal(mu, sd)`, `student_t(df, mu, sd)`,
#' `cauchy(mu, sd)`, `exponential(rate)`, `lkj(eta)`; an empty string
#' sets bounds only.
#'
#' Classes and their scales:
#' - `"b"`: population-level coefficients of `dpar` (default: the
#'   location parameters), excluding the intercept; narrow to one
#'   coefficient with `coef`. Link scale.
#' - `"Intercept"`: the intercept of `dpar`. Link scale.
#' - `"sd"`: random-effect standard deviations (and smoothing SDs), on
#'   the NATURAL sd scale with the log-Jacobian applied, so
#'   `set_prior("exponential(1)", class = "sd")` means what it says;
#'   narrow with `group`.
#' - `"cor"`: the CORRELATION of a random-effect block, as a whole.
#'   `lkj(eta)` only, and it addresses a BLOCK the way class `"sd"`
#'   does, by `group`; `set_prior("lkj(2)", class = "cor")` covers every
#'   correlated block of the model, which is brms's spelling. See
#'   The LKJ prior below.
#' - `"theta"`: raw internal covariance parameters (escape hatch).
#'
#' When priors overlap, later specifications override earlier ones, so
#' put class-wide priors first and coefficient-specific ones after. A
#' class `"theta"` prior on a position an earlier `"cor"` prior covers
#' replaces that whole LKJ term, and the other way round, so "later
#' wins" holds between the two spellings as well.
#' `lb`/`ub` become hard bounds (via Stan's constrained transforms in
#' [frm_sample()]); for class `"sd"` they apply on the sd scale.
#'
#' @section The LKJ prior:
#' `lkj(eta)` is the density `det(C)^(eta - 1)` over a block's
#' correlation matrix `C`, normalized: `eta = 1` is uniform over
#' correlation matrices, larger `eta` concentrates toward the identity.
#' frmtmb holds a correlation as an unconstrained row-normalized
#' Cholesky parameter rather than as `C`, so the density is carried onto
#' those parameters with the exact Jacobian of that map (the derivation
#' is in the source of `R/priors.R`; `tests/testthat/test-lkj.R` checks
#' the sampled correlations against the closed-form LKJ marginals). The
#' prior a FLAT correlation parameter carries instead is
#' `(1 - rho^2)^(-3/2)`, which is improper.
#'
#' It fits `us()` and `gr(cov = )` blocks of two or more terms, which
#' hold a whole correlation matrix, and the one-parameter structures
#' `cs()`, `ar1()` and `hetar1()`, whose single bounded correlation
#' takes the LKJ marginal `(1 - rho^2)^(eta - 1)` with that structure's
#' own Jacobian. A `cs()` correlation is bounded below at `-1/(d - 1)`,
#' where a compound-symmetric matrix stops being positive definite, and
#' the density is renormalized over that window. `toep()` is refused:
#' its parameterization is not positive definite everywhere, so it has
#' no correlation matrix to put a density on.
#'
#' @section Nonlinear parameters:
#' `nlpar` addresses one parameter of an `nl = TRUE` formula, brms's
#' spelling: `set_prior("normal(5000, 1000)", nlpar = "ult")`, or
#' `prior(normal(5000, 1000), nlpar = "ult")`. Class `"b"` there covers
#' EVERY coefficient of that parameter, its intercept included, because
#' a nonlinear parameter's sub-formula is not centered and brms holds
#' its intercept in the same coefficient vector as its slopes. That is
#' why the vignette spelling above lands on `ult_(Intercept)` rather
#' than on nothing. Narrow to one column with `coef` (`"Intercept"` and
#' `"(Intercept)"` both name the intercept), or write
#' `class = "Intercept", nlpar = "ult"`, which is frmtmb's spelling of
#' the same slot. `nlpar` narrows classes `"sd"` and `"cor"` to the
#' random-effect blocks of that parameter as well.
#'
#' A prior with a location places [frm()]'s `start` for a nonlinear
#' parameter. `normal()`, `student_t()` and `cauchy()` all carry one,
#' and where `start` does not set a nonlinear coefficient, that
#' coefficient begins at the prior's location, reported in a message.
#' Other parameters keep their usual starts: a prior is a penalty, not
#' a claim about where to begin. Without a located prior a nonlinear
#' model still needs `start`, because `frm()` evaluates the objective
#' AT the starting values; [par_template()] names them.
#'
#' `resp` picks one response of a multivariate model; the default
#' priors of [frm_sample()] still stay off there (see its Default
#' priors section), so a multivariate model's priors are the ones
#' written by hand.
#'
#' brms's `tag` and `check` have no counterpart: `tag` names a prior
#' for reuse inside a Stan program, and `check` passes an unchecked
#' string through to one. frmtmb compiles no Stan program, so both are
#' omitted rather than accepted and ignored.
#'
#' @param prior Distribution string, e.g. `"normal(0, 5)"`, or a
#'   [prior_normal()]/[prior_t()]/[prior_lkj()] object, or `""` for
#'   bounds only.
#' @param class `"b"`, `"Intercept"`, `"sd"`, `"cor"`, or `"theta"`.
#' @param coef Restrict to one coefficient (classes `"b"`/`"Intercept"`).
#' @param group Restrict class `"sd"` or `"cor"` to one grouping factor.
#' @param resp Response of a multivariate model.
#' @param dpar Distributional parameter (default: the location
#'   parameters).
#' @param nlpar Nonlinear parameter of an `nl = TRUE` formula. See
#'   Nonlinear parameters.
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
#'   `"b"`, `"Intercept"`, `"sd"`, `"cor"` and `"theta"`, so an
#'   unexpected class errors and names the permitted values. The one
#'   distribution that belongs to a single class, `lkj()`, is checked
#'   against it in both directions.
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
#' fit <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, prior = pr)
#' fixef(fit)$mu
#' # the tight prior on z shrinks it toward zero
#' fixef(frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd))$mu
#'
#' # an empty distribution string sets a hard bound only
#' set_prior("", class = "b", coef = "x", lb = 0)
#'
#' # class "cor" addresses a correlated block as a whole, brms's
#' # spelling; eta > 1 pulls the correlation toward zero
#' dd$z <- rnorm(100)
#' dd$y2 <- dd$y + rnorm(10, 0, 0.6)[dd$g] * dd$z
#' fitc <- frm(bf(y2 ~ x + z + (z | g)) + gaussian(), data = dd,
#'             prior = set_prior("lkj(4)", class = "cor"))
#' VarCorr(fitc)
#'
#' # get_prior() shows which rows a design offers
#' get_prior(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
#'
#' # prior() quotes its first argument, brms's spelling, and reaches
#' # the same machinery
#' prior(normal(0, 1), class = "b")
#' @export
set_prior <- function(prior = "", class = "b", coef = "", group = "",
                      resp = "", dpar = "", nlpar = "", lb = NA,
                      ub = NA) {
  dist <- parse_prior_dist(prior)
  # brms carries lb/ub as STRINGS (its prior frame is all character),
  # and prior() deparses everything it is given, so a bound arrives
  # here as "0" as often as 0. Normalizing once means the comparisons
  # downstream are numeric, where `"0" > 0` would have been a string
  # comparison that quietly answered FALSE
  lb <- parse_prior_bound(lb, "lb")
  ub <- parse_prior_bound(ub, "ub")
  if (is.null(dist) && is.na(lb) && is.na(ub)) {
    stop("set_prior() needs a distribution, bounds, or both",
         call. = FALSE)
  }
  class <- match.arg(class, c("b", "Intercept", "sd", "cor", "theta"))
  # lkj is a density over a whole correlation matrix, so it has no
  # meaning on a single coefficient or standard deviation, and class
  # "cor" has no meaning without it: neither mistake can produce a
  # silently different model
  is_lkj <- identical(dist$kind, "lkj")
  if (is_lkj && !identical(class, "cor")) {
    stop("lkj() is a density over a whole correlation matrix; it ",
         "belongs to class = \"cor\" (got class = \"", class, "\")",
         call. = FALSE)
  }
  if (identical(class, "cor") && !is_lkj) {
    stop("class = \"cor\" takes an lkj() prior, e.g. ",
         "set_prior(\"lkj(2)\", class = \"cor\"): it addresses a ",
         "block's whole correlation, which no per-parameter ",
         "distribution describes", call. = FALSE)
  }
  if (identical(class, "cor") && (!is.na(lb) || !is.na(ub))) {
    # a bound belongs to ONE parameter, and this class names a whole
    # correlation matrix; accepting it here would silently drop it
    stop("class = \"cor\" takes no lb/ub: the bound would apply to a ",
         "whole correlation matrix. Bound one parameter at a time with ",
         "class = \"theta\"", call. = FALSE)
  }
  # a nonlinear parameter is addressed by nlpar and a distributional
  # one by dpar; frmtmb's frame gives each its own linear predictor, so
  # naming both at once says one thing twice and could say two
  # different things, which is a question about intent rather than a
  # setting to resolve
  if (nzchar(nlpar) && nzchar(dpar)) {
    stop("set_prior() takes `dpar` or `nlpar`, not both: each ",
         "nonlinear parameter has its own linear predictor here, so ",
         "nlpar = \"", nlpar, "\" already names one slot", call. = FALSE)
  }
  spec <- list(dist = dist, class = class, coef = coef, dpar = dpar,
               group = group, resp = resp, nlpar = nlpar, lb = lb,
               ub = ub)
  structure(list(spec), class = "frmtmb_priorlist")
}

#' A hard bound as one number or `NA`, from the number, the string
#' brms's prior frame stores, or the deparsed constant [prior()]
#' produces.
#'
#' @noRd
parse_prior_bound <- function(x, arg) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  ok <- length(x) == 1L && (is.character(x) || is.numeric(x) ||
                              is.logical(x))
  if (!ok) {
    stop("`", arg, "` must be a single number or NA, not ", arg_desc(x),
         ": a bound here is a hard box constraint on one parameter, ",
         "not a Stan expression", call. = FALSE)
  }
  if (is.na(x)) return(NA_real_)
  if (is.character(x)) {
    if (!nzchar(x)) return(NA_real_)
    v <- suppressWarnings(as.numeric(x))
    if (is.na(v)) {
      stop("`", arg, "` = ", encodeString(x, quote = "\""),
           " is not a number: a bound here is a hard box constraint ",
           "on one parameter, and only a constant can be one",
           call. = FALSE)
    }
    return(v)
  }
  as.numeric(x)
}

#' Set up priors with brms's quoting spelling
#'
#' `prior()` is [set_prior()] with the distribution given UNQUOTED, as
#' brms's `prior()` takes it: `prior(normal(5000, 1000), nlpar = "ult")`
#' is `set_prior("normal(5000, 1000)", nlpar = "ult")`. Every argument
#' is deparsed rather than evaluated, so `class = b` and `class = "b"`
#' mean the same thing, and a variable holding a distribution is
#' deparsed to its NAME rather than its value: use [prior_string()] to
#' build a prior from strings computed at run time.
#'
#' `prior_()` takes one-sided formulas, calls, names or constants
#' (`prior_(~normal(0, 10), class = ~b)`) and `prior_string()` takes
#' plain strings; both exist so that priors can be built
#' programmatically, and both are brms's.
#'
#' A frmtmb prior and a brms prior are different objects, and with brms
#' attached after frmtmb its `prior()` masks this one. Nothing breaks:
#' [frm()] and [frm_sample()] accept a `brmsprior` object and translate
#' its rows, so `c(prior(...), prior(...))` copied out of a brms script
#' works whichever `prior()` was in scope.
#'
#' @inheritParams set_prior
#' @param ... Any of [set_prior()]'s remaining arguments: `class`,
#'   `coef`, `group`, `resp`, `dpar`, `nlpar`, `lb`, `ub`.
#' @return A `frmtmb_priorlist`.
#' @examples
#' # the brms nonlinear vignette's spelling
#' prior(normal(5000, 1000), nlpar = "ult")
#'
#' # combine with c() or `+`, as with set_prior()
#' c(prior(normal(1, 2), nlpar = "omega"),
#'   prior(normal(45, 10), nlpar = "theta"))
#'
#' # the programmatic spellings
#' prior_(~normal(0, 10), class = ~b)
#' prior_string(paste0("normal(0, ", 2 * 5, ")"), class = "b")
#' @export
prior <- function(prior, ...) {
  cl <- as.list(match.call()[-1L])
  do.call(set_prior, lapply(cl, deparse_prior_arg),
          envir = parent.frame())
}

#' @rdname prior
#' @export
prior_ <- function(prior, ...) {
  cl <- c(list(prior = prior), list(...))
  do.call(set_prior, lapply(cl, deparse_prior_value))
}

#' @rdname prior
#' @export
prior_string <- function(prior, ...) set_prior(prior, ...)

#' brms's `deparse_no_string()`: a character argument is already the
#' string `set_prior()` wants, and anything else is the user's
#' unevaluated code, which is deparsed rather than evaluated.
#'
#' @noRd
deparse_prior_arg <- function(x) {
  if (is.character(x)) x else paste(deparse(x), collapse = "")
}

#' The `prior_()` variant, which reads VALUES rather than unevaluated
#' arguments: a one-sided formula gives up its right-hand side, and a
#' call, name or constant is deparsed as it stands.
#'
#' @noRd
deparse_prior_value <- function(x) {
  if (inherits(x, "formula") && length(x) == 2L) {
    return(paste(deparse(x[[2L]]), collapse = ""))
  }
  if (is.character(x)) return(x)
  if (is.call(x) || is.name(x) || is.atomic(x)) {
    return(paste(deparse(x), collapse = ""))
  }
  stop("prior_() takes one-sided formulas, calls, names or constants; ",
       "got ", arg_desc(x), ". prior_string() takes plain strings",
       call. = FALSE)
}

#' Whatever a `prior =` argument turned out to be, as
#' something the resolver understands.
#'
#' A brms `brmsprior` is a data frame of prior/class/coef/group/resp/
#' dpar/nlpar/lb/ub strings, which is exactly `set_prior()`'s vocabulary
#' written down, so it is TRANSLATED rather than refused: with brms
#' attached its `prior()` masks frmtmb's, and a ported script's
#' `c(prior(...), prior(...))` then arrives here as a brms object
#' through no fault of the caller. The legacy named list and the
#' `"flat"` string pass through untouched.
#'
#' @noRd
as_priorlist <- function(x) {
  if (!inherits(x, "brmsprior")) return(x)
  rows <- as.data.frame(x, stringsAsFactors = FALSE)
  chr <- function(nm, i) {
    v <- if (nm %in% names(rows)) rows[[nm]][i] else ""
    if (is.na(v)) "" else as.character(v)
  }
  bnd <- function(nm, i) {
    if (!nm %in% names(rows)) return(NA)
    v <- rows[[nm]][i]
    # brms's frame is all character, and an unset bound is spelled as
    # the empty string in some rows and as NA in others
    if (is.na(v) || (is.character(v) && !nzchar(v))) NA else v
  }
  out <- list()
  dropped <- 0L
  for (i in seq_len(nrow(rows))) {
    dist <- chr("prior", i)
    lb <- bnd("lb", i)
    ub <- bnd("ub", i)
    cls <- chr("class", i)
    # get_prior()/default_prior() rows with an empty `prior` are
    # "this slot exists and is flat", not a prior to apply
    if (!nzchar(dist) && is.na(lb) && is.na(ub)) next
    # a row brms filled in itself is brms's default, and frmtmb chooses
    # its own (see the Default priors section of frm_sample()); keeping
    # both would apply two densities to one parameter
    if (identical(chr("source", i), "default")) {
      dropped <- dropped + 1L
      next
    }
    if (nzchar(chr("tag", i))) {
      stop("A brms prior with tag = \"", chr("tag", i), "\" names a ",
           "prior for reuse inside a Stan program, which frmtmb does ",
           "not build. Drop the tag", call. = FALSE)
    }
    check_brms_prior_class(cls, dist)
    one <- tryCatch(
      set_prior(dist, class = cls, coef = chr("coef", i),
                group = chr("group", i), resp = chr("resp", i),
                dpar = chr("dpar", i), nlpar = chr("nlpar", i),
                lb = lb, ub = ub),
      error = function(e) {
        stop("brms prior row ", i, " (", dist, ", class = \"", cls,
             "\") does not translate: ", conditionMessage(e),
             call. = FALSE)
      })
    out[[length(out) + 1L]] <- unclass(one)[[1L]]
  }
  if (dropped) {
    message("Translating a brms prior: dropped ", dropped,
            " row(s) brms had filled in as its own defaults. frmtmb ",
            "chooses defaults itself on the frm_sample() formula path ",
            "and applies none on the frm() path; write the ones you ",
            "want as prior() rows")
  }
  if (!length(out)) return(NULL)
  structure(out, class = "frmtmb_priorlist")
}

#' The brms classes a translation can honor.
#'
#' brms's class vocabulary is wider than frmtmb's, and the extra names
#' would translate into something DIFFERENT rather than into nothing,
#' which is why each is refused by name rather than passed to
#' `match.arg()`. `theta` is the trap: brms's is a mixture proportion
#' and frmtmb's is the raw internal covariance vector, two unrelated
#' sets of parameters that share the word. A distributional parameter's
#' own class (`sigma`, `shape`, ...) is the near miss: frmtmb's nearest
#' spelling puts that density on the LINK scale, where brms puts it on
#' the parameter itself, so it is a different prior rather than the
#' same one written differently.
#'
#' @noRd
check_brms_prior_class <- function(cls, dist) {
  if (cls %in% c("b", "Intercept", "sd", "cor")) return(invisible(cls))
  hint <- if (identical(cls, "theta")) {
    paste0("brms's \"theta\" is a mixture proportion, and frmtmb's is ",
           "the raw internal covariance vector: the word names two ",
           "unrelated sets of parameters, so the row cannot be carried ",
           "over. ")
  } else {
    paste0("frmtmb's classes are b, Intercept, sd, cor and theta. If ",
           "this names a distributional parameter, its frmtmb spelling ",
           "is class = \"Intercept\", dpar = \"", cls, "\", and that ",
           "density sits on the LINK scale where brms puts it on ", cls,
           " itself. ")
  }
  stop("A brms prior with class = \"", cls, "\" (", dist, ") has no ",
       "faithful frmtmb spelling. ", hint,
       "Write the prior you mean with set_prior() directly",
       call. = FALSE)
}

#' Turn a brms-style prior string such as `"normal(0, 5)"` into a
#' `frmtmb_prior` object, or `NULL` for the empty string. It errors on an
#' unparsable string, on non-numeric arguments, and on the wrong number
#' of arguments for the named distribution.
#'
#' @noRd
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
  pars <- as.numeric(strsplit(m[3], ",", fixed = TRUE)[[1]])
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
    lkj = {
      stopifnot(length(pars) == 1)
      prior_lkj(pars[1])
    },
    stop("Unsupported prior distribution '", kind,
         "' (supported: normal, student_t, cauchy, exponential, lkj)",
         call. = FALSE)
  )
}

#' The location of a prior distribution, or `NA` where it has none.
#' `normal` and `student_t` (which `cauchy` parses into) carry one;
#' `exponential` and `lkj` do not. Read by the nonlinear start
#' placement in R/par-template.R.
#'
#' @noRd
prior_dist_location <- function(dist) {
  loc <- switch(dist$kind %||% "", normal = , t = dist$location, NULL)
  if (is.null(loc) || length(loc) != 1L || !is.finite(loc)) {
    return(NA_real_)
  }
  as.numeric(loc)
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
    d <- if (is.null(s$dist)) "(bounds only)" else {
      # brms spelling on the way out as well as on the way in, so a
      # printed prior can be pasted back into set_prior()
      kind <- if (identical(s$dist$kind, "t")) "student_t" else s$dist$kind
      paste0(kind, "(", paste(unlist(s$dist[-1]), collapse = ", "), ")")
    }
    cat(d, " class=", s$class,
        if (nzchar(s$coef)) paste0(" coef=", s$coef),
        if (nzchar(s$dpar)) paste0(" dpar=", s$dpar),
        if (nzchar(s$nlpar %||% "")) paste0(" nlpar=", s$nlpar),
        if (nzchar(s$resp %||% "")) paste0(" resp=", s$resp),
        if (nzchar(s$group)) paste0(" group=", s$group),
        if (isTRUE(s$natural)) " scale=natural",
        if (!is.na(s$lb)) paste0(" lb=", s$lb),
        if (!is.na(s$ub)) paste0(" ub=", s$ub), "\n", sep = "")
  }
  ov <- attr(x, "overrides")
  if (length(ov)) {
    cat("plus internal-scale overrides on: ",
        paste(names(ov), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Enumerate the targetable prior slots
#'
#' The [set_prior()] counterpart of brms's `get_prior()`: one row per
#' slot a prior can target, with the class/coef/dpar/group values to
#' pass to `set_prior()`. The default in every slot is flat (this is
#' maximum likelihood until priors are set; the formula route of
#' [frm_sample()] has its own brms defaults, which
#' `prior_summary()` reports). Classes `"sd"` and `"cor"` are targeted
#' by `group` and `nlpar`; class `"theta"` rows name the raw internal
#' covariance parameters (escape hatch, including correlations one at a
#' time).
#'
#' A nonlinear parameter's coefficients are listed under class `"b"`
#' with its name in the `nlpar` column, the intercept among them, which
#' is how brms lists them and what [set_prior()] addresses (see its
#' Nonlinear parameters section).
#'
#' @param formula A `bf()` formula (with family), a plain formula, or
#'   an already fitted `frmtmb_fit`.
#' @param data A data frame of model data (ignored when `formula` is a
#'   fit).
#' @param family Family, when `formula` does not carry one.
#' @param data2 Structural objects, as in [frm()] (ignored when
#'   `formula` is a fit, which carries its own).
#' @return A data frame with columns `prior`, `class`, `coef`,
#'   `group`, `dpar`, `nlpar`, `resp`, `lb`, `ub`.
#' @examples
#' dd <- data.frame(y = rnorm(60), x = rnorm(60),
#'                  g = factor(rep(1:6, 10)))
#' get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' @export
get_prior <- function(formula, data = NULL, family = NULL,
                      data2 = list()) {
  if (inherits(formula, "frmtmb_fit")) {
    spec <- formula$spec
    frame <- formula$frame
  } else {
    bform <- resolve_deferred_families(as_bform(formula, family), data)
    spec <- parse_spec(bform)
    frame <- assemble_frame(spec, data, data2 = data2)
  }

  multi <- length(spec$responses) > 1L
  rows <- list()
  add <- function(class, coef = "", group = "", dpar = "", nlpar = "",
                  resp = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      prior = "(flat)", class = class, coef = coef, group = group,
      dpar = dpar, nlpar = nlpar, resp = resp, lb = NA_real_,
      ub = NA_real_
    )
  }

  for (lp in frame$linpreds) {
    if (!is.null(lp$constant) || !is.null(lp$nl_body)) next
    rspec <- spec$responses[[lp$resp]]
    # location dpars are the default target (dpar = ""), matching
    # set_prior()'s resolution
    dpar_lab <- if (lp$dpar %in% rspec$primary_dpars) "" else lp$dpar
    resp_lab <- if (multi) lp$resp else ""
    nl_lab <- if (lp$dpar %in% (rspec$nlpars %||% character(0))) {
      lp$dpar
    } else {
      ""
    }
    cn <- colnames(lp$X)
    # a nonlinear parameter's sub-formula is not centered, so its whole
    # coefficient vector including the intercept is class "b" - brms's
    # own listing, and what set_prior(nlpar =) addresses
    if (nzchar(nl_lab)) {
      if (length(cn)) {
        add("b", dpar = dpar_lab, nlpar = nl_lab, resp = resp_lab)
        for (co in cn) {
          add("b", coef = co, dpar = dpar_lab, nlpar = nl_lab,
              resp = resp_lab)
        }
      }
      next
    }
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

  sd_rows <- list()
  cor_rows <- list()
  for (bk in frame$re_blocks) {
    key <- list(group = bk$group_name, nlpar = block_nlpar(spec, frame, bk),
                resp = if (multi) block_resp(frame, bk) else "")
    sd_i <- covstruct_registry[[bk$covstruct]]$sd_idx(bk$dim)
    if (length(sd_i)) sd_rows[[length(sd_rows) + 1L]] <- key
    if (!bk$covstruct %in% names(lkj_refusals) &&
        !is.null(block_cor_spec(bk))) {
      cor_rows[[length(cor_rows) + 1L]] <- key
    }
  }
  for (cl in c("sd", "cor")) {
    ks <- if (identical(cl, "sd")) sd_rows else cor_rows
    if (!length(ks)) next
    add(cl)
    for (k in ks) {
      add(cl, group = k$group, nlpar = k$nlpar, resp = k$resp)
    }
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

#' The linear predictors a random-effect block draws its columns from.
#' A block merges the components that share an `|ID|` key, and those can
#' come from different predictors, so this is a set rather than one
#' value; `bk$dpar` is the first component's and is the fallback for a
#' block whose components predate the key.
#'
#' @noRd
block_linpreds <- function(frame, bk) {
  keys <- vapply(bk$components %||% list(),
                 function(cp) cp$lp_key %||% "", "")
  lps <- frame$linpreds[keys[nzchar(keys)]]
  Filter(Negate(is.null), lps)
}

#' The nonlinear parameter a block belongs to, or `""`. Blocks that
#' straddle several are left unlabeled rather than assigned to one.
#'
#' @noRd
block_nlpar <- function(spec, frame, bk) {
  lps <- block_linpreds(frame, bk)
  np <- unique(vapply(lps, function(lp) {
    nl <- spec$responses[[lp$resp]]$nlpars %||% character(0)
    if (lp$dpar %in% nl) lp$dpar else ""
  }, ""))
  np <- np[nzchar(np)]
  if (length(np) == 1L) np else ""
}

#' The response a block belongs to, or `""` when it straddles several.
#'
#' @noRd
block_resp <- function(frame, bk) {
  rs <- unique(vapply(block_linpreds(frame, bk), `[[`, "", "resp"))
  if (length(rs) == 1L) rs else ""
}

#' Whether a class `"sd"` / `"cor"` specification addresses this block.
#' `group` has always narrowed here; `resp`, `dpar` and `nlpar` narrow
#' the same way, which is what lets a nonlinear model's per-parameter
#' variance components be priored one at a time.
#'
#' @noRd
block_addressed <- function(spec, frame, bk, s) {
  if (nzchar(s$group) && !identical(bk$group_name, s$group)) return(FALSE)
  want_np <- nzchar(s$nlpar %||% "")
  want_rs <- nzchar(s$resp %||% "")
  want_dp <- nzchar(s$dpar)
  if (!want_np && !want_rs && !want_dp) return(TRUE)
  lps <- block_linpreds(frame, bk)
  if (want_rs && !s$resp %in% vapply(lps, `[[`, "", "resp")) return(FALSE)
  dp <- vapply(lps, `[[`, "", "dpar")
  if (want_np && !identical(block_nlpar(spec, frame, bk), s$nlpar)) {
    return(FALSE)
  }
  if (want_dp && !s$dpar %in% dp) return(FALSE)
  TRUE
}

#' Where a specification was pointing, for the refusals that report a
#' target no design offers. Written once so both refusals name the same
#' fields in the same order.
#'
#' @noRd
spec_target <- function(s) {
  paste0("class=", s$class,
         if (nzchar(s$coef)) paste0(", coef=", s$coef),
         if (nzchar(s$group)) paste0(", group=", s$group),
         if (nzchar(s$resp %||% "")) paste0(", resp=", s$resp),
         if (nzchar(s$dpar)) paste0(", dpar=", s$dpar),
         if (nzchar(s$nlpar %||% "")) paste0(", nlpar=", s$nlpar))
}

#' The nonlinear parameters a model declares, for the refusal that has
#' to say what `nlpar` could have named.
#'
#' @noRd
model_nlpars <- function(spec) {
  unique(unlist(lapply(spec$responses, function(r) {
    r$nlpars %||% character(0)
  }))) %||% character(0)
}

#' Copy the bounds of a prior specification onto an existing entry. This
#' is how a later bounds-only specification tightens an entry that an
#' earlier distribution created, instead of being lost.
#'
#' @noRd
entry_bounds <- function(entry, s) {
  if (!is.na(s$lb)) entry$lb <- s$lb
  if (!is.na(s$ub)) entry$ub <- s$ub
  entry
}

#' Resolve a priorlist against a fit: per-parameter prior entries (later
#' specifications override earlier ones) plus named bound vectors. An
#' entry holds `comp`, a scalar `idx`, `dist`, `scale` ("internal" or
#' "sd"), and `lb`/`ub` on the entry's own scale. `frm_simulate()`
#' rejects draws outside those bounds; `frm_sample()` uses the named
#' bound vectors instead.
#'
#' @noRd
resolve_priorlist <- function(fit, pl) {
  frame <- fit$frame
  assigned <- list()   # key "comp.idx" -> entry
  lower <- c()
  upper <- c()
  # a class "cor" entry covers SEVERAL theta positions at once, so the
  # key names them all; a one-position entry keys exactly as before
  nm_of <- function(comp, idx) paste0(comp, ".", paste(idx, collapse = ","))

  # Assigning over a position that a joint (class "cor") entry already
  # covers has to RETIRE that entry, and so does a joint entry covering
  # positions that per-parameter entries claimed. Otherwise both
  # densities would be added. With that, "later wins" reads the same
  # across the two spellings.
  claim <- function(comp, idx) {
    drop <- vapply(assigned, function(e) {
      identical(e$comp, comp) && length(intersect(e$idx, idx)) > 0L
    }, TRUE)
    assigned <<- assigned[!drop]
  }

  nlpars <- model_nlpars(fit$spec)

  target_coefs <- function(s) {
    # (comp, idx, name) triplets for classes b / Intercept
    out <- list()
    want_np <- nzchar(s$nlpar %||% "")
    if (want_np && !s$nlpar %in% nlpars) {
      stop("nlpar = \"", s$nlpar, "\" names no nonlinear parameter of ",
           "this model. It has ",
           if (length(nlpars)) paste(nlpars, collapse = ", ") else
             "none (write nl = TRUE in bf() to declare them)",
           ". A distributional parameter is addressed with dpar =",
           call. = FALSE)
    }
    for (lp in frame$linpreds) {
      if (!is.null(lp$constant) || !is.null(lp$nl_body)) next
      rspec <- fit$spec$responses[[lp$resp]]
      if (nzchar(s$resp %||% "") && !identical(lp$resp, s$resp)) next
      is_loc <- lp$dpar %in% rspec$primary_dpars
      if (want_np) {
        if (!identical(lp$dpar, s$nlpar)) next
      } else if (nzchar(s$dpar)) {
        if (!identical(lp$dpar, s$dpar)) next
      } else if (!is_loc) next
      cn <- colnames(lp$X)
      pick <- if (s$class == "Intercept") {
        which(cn == "(Intercept)")
      } else if (nzchar(s$coef)) {
        # brms writes an intercept as "Intercept"; the design matrix
        # spells it "(Intercept)", and both name the same column
        which(cn == s$coef | par_name_bare(cn) == par_name_bare(s$coef))
      } else if (want_np) {
        # a nonlinear parameter's sub-formula is NOT centered, so its
        # intercept sits in the same coefficient vector as its slopes
        # and class "b" covers it. This is the whole reason brms's
        # prior(normal(5000, 1000), nlpar = "ult") lands on an
        # intercept-only nonlinear parameter
        seq_along(cn)
      } else {
        which(cn != "(Intercept)")
      }
      for (k in pick) {
        out[[length(out) + 1L]] <- list(comp = lp$par,
                                        idx = lp$idx[k], name = cn[k])
      }
    }
    if (!length(out) &&
          (nzchar(s$coef) || s$class == "Intercept" || want_np ||
             nzchar(s$resp %||% ""))) {
      stop("Prior target not found (", spec_target(s), ")",
           call. = FALSE)
    }
    out
  }

  for (s in unclass(pl)) {
    if (s$class %in% c("b", "Intercept")) {
      # `natural` puts the prior on exp(coefficient) with the same
      # log-Jacobian class "sd" uses, which is what a log-linked
      # dispersion intercept needs to carry brms's half-t on sigma
      # itself rather than on log sigma. Only the default-prior builder
      # of frm_sample() sets it; a set_prior() spec never has the field
      # and reads as internal, exactly as before.
      sc <- if (isTRUE(s$natural)) "sd" else "internal"
      for (tg in target_coefs(s)) {
        key <- nm_of(tg$comp, tg$idx)
        if (!is.null(s$dist)) {
          assigned[[key]] <- list(comp = tg$comp, idx = tg$idx,
                                  dist = s$dist, scale = sc,
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
        if (!block_addressed(fit$spec, frame, bk, s)) next
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
        stop("No random-effect SDs match ", spec_target(s),
             call. = FALSE)
      }
    } else if (s$class == "cor") {
      hit <- FALSE
      refused <- character(0)
      for (bk in frame$re_blocks) {
        if (!block_addressed(fit$spec, frame, bk, s)) next
        cs <- bk$covstruct
        if (cs %in% names(lkj_refusals)) {
          refused <- c(refused, paste0(bk$term_label, " [", cs, "]: ",
                                       unname(lkj_refusals[[cs]])))
          next
        }
        spec <- block_cor_spec(bk)
        if (is.null(spec)) next
        hit <- TRUE
        idx <- bk$theta_idx[spec$idx]
        claim("theta", idx)
        assigned[[nm_of("theta", idx)]] <-
          list(comp = "theta", idx = idx,
               dist = lkj_dist(s$dist$eta, spec), scale = "internal",
               lb = NA, ub = NA)
      }
      if (!hit) {
        # a refused structure is named with its reason; otherwise the
        # model simply has no correlation to prior, and the message says
        # what it does have
        have <- unique(vapply(frame$re_blocks, function(bk) {
          paste0(bk$term_label, " [", bk$covstruct, "]")
        }, ""))
        stop("No random-effect correlations match ", spec_target(s),
             ". ",
             if (length(refused)) {
               paste0("No LKJ density fits ",
                      paste(refused, collapse = "; "))
             } else if (length(have)) {
               paste0("These blocks have no correlation parameter: ",
                      paste(have, collapse = ", "))
             } else {
               "This model has no random-effect blocks"
             }, call. = FALSE)
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
          claim("theta", i)
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

#' Log density of one prior entry value (AD-safe), with the sd-scale
#' change of variables where requested.
#'
#' @noRd
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
    exponential = log(dist$rate) - dist$rate * x,
    # a JOINT density over a whole correlation, so `x` is the block's
    # correlation segment and the value is one number, not one per
    # element (see lkj_logdens)
    lkj = lkj_logdens(x, dist)
  )
  ld + jac
}

#' @param eta LKJ shape. `1` is uniform over correlation matrices,
#'   larger values concentrate toward the identity, and `0 < eta < 1`
#'   pushes toward the boundary.
#' @rdname frmtmb-priors
#' @export
prior_lkj <- function(eta = 1) {
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) ||
      eta <= 0) {
    stop("prior_lkj(eta =) takes one finite positive number; eta = 1 ",
         "is uniform over correlation matrices", call. = FALSE)
  }
  structure(list(kind = "lkj", eta = eta), class = "frmtmb_prior")
}

# ------------------------------------------------- the LKJ prior ------
#
# THE DENSITY, ON FRMTMB'S OWN PARAMETERS. LKJ(eta) is the density
# `p(C) = c_d(eta)^-1 det(C)^(eta - 1)` over d x d correlation matrices;
# `eta = 1` is uniform over them. frmtmb never holds `C`. It holds `t`,
# the strictly-lower entries of a unit-diagonal lower-triangular matrix
# whose rows are then normalized (us_chol_L()), so the prior has to be
# carried onto `t` with the Jacobian of that map, and the result is what
# is implemented here. The derivation, in three steps:
#
# 1. ON THE CHOLESKY FACTOR. With `C = L L'`, `L` lower-triangular with
#    positive diagonal and unit-norm rows, the Jacobian of the map from
#    the strict lower triangle of `C` to that of `L` is
#    `prod_{i>=2} L_ii^(d - i)`, so
#      p(L) = c_d(eta)^-1 prod_{i>=2} L_ii^(d - i + 2 eta - 2),
#    using `det(C) = prod L_ii^2`. That is Stan's
#    `lkj_corr_cholesky_lpdf`, and it is the standard route.
#
# 2. FROM L TO t, ROW BY ROW. Row `i` of frmtmb's unnormalized matrix is
#    `(t_i, 1, 0, ...)` with `t_i` of length `m = i - 1`, so row `i` of
#    `L` is that vector over `sqrt(1 + ||t_i||^2)` and
#      L_ii = (1 + ||t_i||^2)^(-1/2).
#    The free entries of row `i` of `L` are `u_i = t_i / sqrt(1 + s)`,
#    `s = ||t_i||^2`, whose Jacobian matrix is
#    `(1 + s)^(-1/2) (I - t_i t_i' / (1 + s))`, with determinant
#      (1 + s)^(-m/2) * (1 - s/(1 + s)) = (1 + s)^(-(m + 2)/2)
#                                       = L_ii^(m + 2) = L_ii^(i + 1).
#
# 3. THE PRODUCT. Row `i` therefore contributes
#    `L_ii^(d - i + 2 eta - 2) * L_ii^(i + 1) = L_ii^(2 eta + d - 1)`:
#    the exponent is the SAME for every row, and with
#    `log L_ii = -log(1 + ||t_i||^2)/2` the whole log density is
#
#      log p(t) = -(eta + (d - 1)/2) * sum_i log(1 + ||t_i||^2)
#                 - sum_i log Z_i.
#
#    `Z_i` normalizes row `i` on its own, because the rows are
#    independent under LKJ: `p(u_i) ∝ (1 - ||u_i||^2)^(a_i - 1)` on the
#    unit ball of R^m with `a_i = eta + (d - i)/2`, and
#    `int (1 - ||u||^2)^(a-1) du = pi^(m/2) Gamma(a) / Gamma(a + m/2)`,
#    which gives `Z_i = pi^((i-1)/2) Gamma(eta + (d-i)/2) /
#    Gamma(eta + (d-1)/2)`. Their product IS the published `c_d(eta)`
#    (checked to 1e-15 against the LKJ 2009 closed form for d = 2..5,
#    tests/testthat/test-lkj.R).
#
# THE d = 2 CHECK, closed form: `rho = t / sqrt(1 + t^2)`, so
# `p(t) = (1 + t^2)^-(eta + 1/2) / (2^(2 eta - 1) B(eta, eta))`, which at
# `eta = 1` is `(1 + t^2)^(-3/2)` and is uniform on rho. Flat on `t`, by
# the same change of variables, is `(1 - rho^2)^(-3/2)`: improper, all
# its mass at |rho| = 1, which is what the LKJ prior replaces.
#
# THE ONE-PARAMETER STRUCTURES (cs, homcs, ar1, hetar1) hold a single
# bounded correlation instead of a whole matrix, so there is no matrix
# for the density above to be about. They take the d = 2 form,
# `p(rho) ∝ (1 - rho^2)^(eta - 1)`, which is the LKJ marginal, with each
# structure's own Jacobian: `rho = t/sqrt(1 + t^2)` for ar1 (identical to
# the d = 2 case above), and the scaled logistic onto `(-1/(d-1), 1)` for
# cs, whose normalizer picks up the mass LKJ puts below `-1/(d-1)`.

#' `sum_i log Z_i`: the log of the LKJ normalizing constant for
#' dimension `d`, assembled from the per-row constants of the derivation
#' above.
#'
#' @noRd
lkj_lognorm <- function(eta, d) {
  if (d < 2L) return(0)
  i <- seq_len(d - 1L) + 1L
  sum((i - 1) / 2 * log(pi) + lgamma(eta + (d - i) / 2) -
        lgamma(eta + (d - 1) / 2))
}

#' The positions of the correlation segment belonging to each ROW of the
#' Cholesky factor. `L[lower.tri(L)] <- t` fills column-major, so the
#' entries of one row are not contiguous.
#'
#' @noRd
lkj_rows <- function(d) {
  ii <- row(diag(d))[lower.tri(diag(d))]
  lapply(seq_len(d - 1L) + 1L, function(i) which(ii == i))
}

#' The internal prior object the objective evaluates: the user's `eta`
#' plus everything about the block's map that does not depend on the
#' parameters, computed once here so the taped density is arithmetic
#' only.
#'
#' @noRd
lkj_dist <- function(eta, spec) {
  d <- spec$d
  out <- list(kind = "lkj", eta = eta, map = spec$kind, d = d)
  if (identical(spec$kind, "chol")) {
    out$rows <- lkj_rows(d)
    out$pow <- eta + (d - 1) / 2
    out$lognorm <- lkj_lognorm(eta, d)
  } else if (identical(spec$kind, "ar1")) {
    out$pow <- eta + 0.5
    out$lognorm <- lkj_lognorm(eta, 2L)
  } else {
    a <- 1 / (d - 1)
    out$a <- a
    # the marginal restricted to (-a, 1), renormalized over that window
    out$lognorm <- lkj_lognorm(eta, 2L) +
      log(1 - stats::pbeta((1 - a) / 2, eta, eta))
  }
  structure(out, class = "frmtmb_prior")
}

#' The LKJ log density at one block's correlation parameters (AD-safe).
#' One number for the whole block, whatever its map.
#'
#' @noRd
lkj_logdens <- function(t, dist) {
  if (identical(dist$map, "chol")) {
    q <- 0
    for (r in dist$rows) q <- q + log(1 + sum(t[r] * t[r]))
    return(-dist$pow * q - dist$lognorm)
  }
  if (identical(dist$map, "ar1")) {
    return(-dist$pow * log(1 + t[1] * t[1]) - dist$lognorm)
  }
  a <- dist$a
  rho <- -a + (1 + a) / (1 + exp(-t[1]))
  # d rho / d t = (rho + a)(1 - rho)/(1 + a), written in rho so that the
  # logistic is computed once
  ld <- log(rho + a) + log(1 - rho) - log(1 + a) - dist$lognorm
  # eta = 1 has no density factor at all, and writing the term anyway
  # would evaluate 0 * log(0) = NaN where the logistic saturates. `eta`
  # is a constant, so this branch is resolved when the tape is built.
  if (dist$eta != 1) ld <- ld + (dist$eta - 1) * log(1 - rho * rho)
  ld
}

#' Accepts the legacy named list of prior objects OR a priorlist; returns
#' entries plus bounds.
#'
#' @noRd
resolve_prior_input <- function(fit, prior) {
  # the argument boundaries translate a brms prior object already; this
  # covers the internal callers that reach the resolver directly
  prior <- as_priorlist(prior)
  if (inherits(prior, "frmtmb_priorlist")) {
    return(resolve_priorlist(fit, prior))
  }
  legacy <- resolve_priors(fit, prior)
  entries <- list()
  for (e in legacy) {
    for (i in e$idx) {
      entries[[length(entries) + 1L]] <-
        list(comp = e$comp, idx = i, dist = e$prior, scale = "internal")
    }
  }
  list(entries = entries, lower = c(), upper = c())
}