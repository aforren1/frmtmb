#' The refusal a function whose `...` would otherwise swallow the
#' retired `priors =` spelling raises instead.
#'
#' The argument is `prior`, brms's spelling and the only one. R's
#' partial matching cannot bind `priors` to it (a longer name is not a
#' prefix of a shorter one), so every direct call fails loudly on its
#' own. The exception is a function whose `...` stands ready to absorb
#' the name in silence and fit the model with no priors at all, and
#' frm_sample() is the only one in either package.
#'
#' @noRd
refuse_retired_priors <- function(dots, what) {
  if ("priors" %in% names(dots)) {
    frm_stop(what, " takes `prior`, not `priors`: the argument follows ",
             "brms's spelling, and this function's `...` would otherwise ",
             "pass the old name through and fit with no priors at all. ",
             "Rename it to `prior`", call. = FALSE)
  }
  # `lower`/`upper` retired in 0.49. Here they would reach tmbstan as
  # unknown sampler options rather than failing as unused arguments, so
  # the model would sample UNBOUNDED in silence
  bad <- intersect(c("lower", "upper"), names(dots))
  if (length(bad)) {
    frm_stop(what, " has no `", paste(bad, collapse = "`/`"),
             "`: a hard bound is written as a prior, ",
             "set_prior(\"\", nlpar = \"la\", lb = 0), which reaches the ",
             "same constrained transform. Every outer parameter has a ",
             "class, down to one internal covariance parameter ",
             "(class = \"theta\", coef = \"theta_1\"). Left in `...` ",
             "these would have gone to the sampler as unknown options and ",
             "the model would have sampled unbounded", call. = FALSE)
  }
  invisible(NULL)
}

# --- the two `control` lists -------------------------------------------
#
# `frm_sample(control =)` is the SAMPLER's control list, which is what
# brms means by the name, so a brms call ports across unchanged. The
# fit-time options frm() takes are `fit_control =`. The two names used
# to be the other way round, and adapt_delta then had no route at all;
# vignette("brms-posterior") carried a row saying so.
#
# Nothing bridges the two names. A frmtmb_control() list arriving in
# `control` is refused rather than reinterpreted, because passing its
# fields on would ask rstan to adapt on `grad_tol`.

#' rstan's own control vocabulary, copied from `rstan:::config_argss()`.
#'
#' Checked here rather than left to rstan because rstan does not raise
#' on an unknown option: it prints the complaint, declines to sample,
#' and returns an EMPTY stanfit. That lands on the "returned no draws"
#' refusal further down, whose text is about a tape calling an external
#' solver, so a typo in `control` would be reported as a solver failure.
#'
#' @noRd
stan_control_names <- c(
  "adapt_engaged", "adapt_gamma", "adapt_delta", "adapt_kappa",
  "adapt_t0", "adapt_init_buffer", "adapt_term_buffer", "adapt_window",
  "stepsize", "stepsize_jitter", "metric", "int_time", "max_treedepth",
  "epsilon", "error"
)

#' Refuse a `frmtmb_control()` list where the sampler's control belongs.
#'
#' frmtmb_control() returns a PLAIN list, so there is no class to
#' dispatch on and the two lists have to be told apart by their field
#' names. That works because the two vocabularies are disjoint: rstan's
#' fifteen names above share none of frmtmb_control()'s fourteen. The
#' refusal names the fields it keyed on, so it reports what it saw
#' rather than asserting what the caller meant.
#'
#' @noRd
check_stan_control <- function(control, what) {
  # an empty list is a legitimate no-op, and it is what building the
  # list programmatically gives when no option was added
  if (is.null(control) || (is.list(control) && !length(control))) {
    return(invisible(NULL))
  }
  if (!is.list(control)) {
    frm_stop(what, " takes `control` as a NAMED list of rstan sampler ",
             "options, control = list(adapt_delta = 0.99), which is the ",
             "spelling brms uses; got an object of class ",
             paste(class(control), collapse = "/"),
             ". The fit-time options of frmtmb_control() are ",
             "`fit_control` now", call. = FALSE)
  }
  nms <- names(control) %||% rep("", length(control))
  if (!all(nzchar(nms))) {
    frm_stop(what, ": every element of `control` names a sampler option, ",
             "and ", sum(!nzchar(nms)), " of ", length(control),
             " have no name. Write control = list(adapt_delta = 0.99); an ",
             "unnamed element cannot reach rstan as anything",
             call. = FALSE)
  }
  fitside <- intersect(nms, names(frmtmb_control()))
  if (length(fitside)) {
    frm_stop(what, ": `control` is the SAMPLER's control list here, as in ",
             "brms (control = list(adapt_delta = 0.99)), and this one ",
             "carries the frmtmb_control() field(s) ",
             paste(fitside, collapse = ", "),
             ". Those are fit-time options and go to `fit_control` now: ",
             "fit_control = frmtmb_control(...). They are not passed on to ",
             "rstan, which would otherwise be asked to adapt on them",
             call. = FALSE)
  }
  bad <- setdiff(nms, stan_control_names)
  if (length(bad)) {
    frm_stop(what, ": rstan has no sampler option named ",
             paste(bad, collapse = ", "), ". The options are ",
             paste(stan_control_names, collapse = ", "),
             ". rstan does not raise on an unknown one; it declines to ",
             "sample and returns an empty fit, so this is refused here ",
             "instead", call. = FALSE)
  }
  invisible(NULL)
}

#' Refuse an abbreviation of `control` that would reach rstan anyway.
#'
#' `control` follows `...` in this function's formals, so R binds only
#' its exact spelling and an abbreviation lands in `...` instead. That
#' is where the reasoning used to stop, and it stopped one step early:
#' `...` is handed to `tmbstan::tmbstan()` and then to
#' `rstan::sampling()`, whose `control` argument comes BEFORE that
#' method's own `...`, so R partial-matches the abbreviation onto it
#' there. Measured: `frm_sample(fit, contro = list(adapt_delta = 0.97))`
#' reached the sampler with adapt_delta 0.97 and never passed
#' `check_stan_control()`, and `contro = frmtmb_control()` reached it
#' too and came back as the "returned no draws ... external solver"
#' report, which is the exact mis-report that check exists to prevent.
#'
#' rstan already rejects a dots name that is not a prefix of one of its
#' own formals ("passing unknown arguments: fit_contro"), so the hole
#' is only the partial-match set, and `control` is the one member of it
#' this package validates itself.
#'
#' @noRd
refuse_partial_control <- function(dots, what) {
  nms <- names(dots)
  if (is.null(nms)) return(invisible(NULL))
  # duplicates.ok, or a second abbreviation in the same call matches
  # nothing and slips through
  hit <- !is.na(pmatch(nms, "control", duplicates.ok = TRUE))
  bad <- nms[nzchar(nms) & hit & nms != "control"]
  if (length(bad)) {
    frm_stop(what, ": `", paste(bad, collapse = "`, `"),
             "` is an abbreviation of `control`, and an abbreviation is ",
             "not harmless here. It lands in `...`, and `...` reaches ",
             "rstan::sampling(), where `control` comes before that ",
             "function's own `...` and partial matching binds it, so the ",
             "list would reach the sampler without being checked. Spell ",
             "it `control`; the fit-time options are `fit_control`",
             call. = FALSE)
  }
  invisible(NULL)
}

# --- the pre-flight against the compatibility registry -----------------
#
# A pair this sampler cannot run is declared once, in the registry, by
# whichever package makes the pair possible, and read here. The
# alternative is a guard written into this file that names a package
# this one does not depend on, and a registry row saying the same thing
# that nothing consults; the row would then be documentation of a
# refusal rather than the refusal itself.
#
# WHAT IT CAN SEE, AND WHY THE ANSWER IS AS NARROW AS IT IS.
#
# What decides whether a model uses a registry feature is the POSITION
# of a call in the formula: `se` inside the bar on the response's left
# is the `se()` addition term, and `se` in a nonlinear body is whatever
# function the user has of that name. The parser owns position, and
# most of what it knows is not carried out of core, so this pre-flight
# decides three ways.
#
#   1. kind "family": matched against the response families by name. A
#      family is named, never called, so no call can be confused with
#      one. It does NOT decompose mixture(): a mixture's family name is
#      the whole string "mixture(gaussian, gaussian)", so a refused
#      `gaussian` row misses a mixture of gaussians. A miss.
#   2. kind "aterm", spelled with parentheses: matched BY POSITION.
#      An addition term is written inside the bar on the left-hand side
#      of the response formula and nowhere else, and the bf() object
#      carries that, so this needs nothing from core and is sound. It
#      is what separates `bf(y | se(v) ~ x)` from
#      `bf(y ~ a * se(x), a ~ 1, nl = TRUE)`, which frmtmb fits.
#   3. any other display name ending in "()": matched BY NAME, on the
#      display name minus its parentheses, over calls OUTSIDE aterm
#      position, and only when a formula actually writes that name and
#      it is not a function of a base-priority R package. Both tests
#      are read off the registry and the installed base packages at
#      call time, not assumed.
#
# WHAT (3) BUYS, EXACTLY, because the first version of this comment
# overstated it. The two tests remove two ENUMERABLE classes of name
# collision: a name two features share, and a name base R also uses.
# They do not touch the class that motivates the whole design, which is
# a call written where the feature cannot occur. `mo`, `ma`, `cosy` and
# `mm` all pass both tests and would all fire on a user's own function
# of that name in a nonlinear body. So (3) is the largest set
# justifiable FROM NAMES ALONE, not a sound one. It is used only where
# position is unavailable, and every remaining hole is a false ALARM
# that would need someone to register a refusal on one of those rows
# first; the only refused row today is `frm_ode()`.
#
# WHY THE DISPLAY NAME AND NOT THE REGISTRY KEY. The key is the parser
# name, which reads better, and the first version used it. It is also
# shared: `mi()`/`mi_pred()`, `gp`/`gp_pred()` and `cs`/`cs_pred()` each
# share one, a collision R/compat.R records as deliberate. Keying on it
# therefore had to degrade to a warning whenever a key was shared, and
# that made a REGISTERED REFUSAL REVOCABLE BY A THIRD PARTY: registering
# a bare `frm_ode` display name, which core accepts, downgraded
# frmtmb.ode's refusal to a warning, and the model then sampled, which
# is the thing that leaves rstan unusable for the rest of the session.
#
# The display name cannot be attacked that way, and the reason is a
# proof rather than a measurement. Two call-shaped rows cannot want the
# same call name, because core refuses a second registration of a
# display name ("One display name carries one KIND", and the kind is as
# load-bearing here as the name, since the route is chosen from both),
# and because dropping a fixed "()" suffix is injective over the names
# that carry it. A bare name cannot compete for a call either, and that
# is a property of the CODE and not of the vocabulary: `ar1` and
# `gr_cov` are bare rows that a formula DOES write as calls, and what
# stops them is that the parentheses test above returns before a bare
# name can reach this route at all.
#
# What matching the display name costs is the three `_pred` specials,
# whose display names no formula writes. They do not silently miss:
# case 3 compares the display name against the registry KEY, which is
# the parser name, and routes the three rows where the two differ to
# the warning. Silence there would be the failing-open shape this file
# exists to close. The filed seam is what would decide them.
#
# Everything else warns and samples: a bare display name, a fitting
# mode, another post-fit method. Warning rather than passing in silence
# keeps a registered refusal from looking like a passed check; warning
# rather than refusing keeps the check from firing on a correct model.
# `frm_features_used(spec)` is filed in dev/extension-gaps-plan.md and
# would close the `special`, `autocor` and `grammar` kinds; case 2
# above is the `aterm` kind closed without it.

#' The calls a model's formulas make, split by whether they sit in
#' ADDITION-TERM position.
#'
#' Returns `list(aterm =, other =)`. An addition term is written inside
#' the bar on the left-hand side of a response formula, `y | se(v) ~ x`,
#' and nowhere else; everything else a formula calls goes in `other`.
#' That split is the one position test available without reading
#' anything from core, and it is what separates the `se()` addition
#' term from a user's own `se()` in a nonlinear body.
#'
#' Walks the `bf()` object, which is where a call-shaped feature is
#' written: the response formula, the dpar and nonlinear formulas of
#' `pforms`, and the `nlf()` bodies. Reached on both of
#' `frm_sample()`'s routes, so the refusal lands before the formula
#' route assembles anything.
#'
#' Only the RESPONSE formula has an aterm side. A `pforms` entry is a
#' distributional or nonlinear parameter's formula and a `|` there is
#' something else, so those contribute to `other` in full.
#'
#' @noRd
formula_calls <- function(x) {
  out <- list(aterm = character(0), other = character(0))
  add <- function(a, b) {
    list(aterm = unique(c(a$aterm, b$aterm)),
         other = unique(c(a$other, b$other)))
  }
  if (inherits(x, "frmtmb_mvformula")) {
    for (f in x[["forms"]]) out <- add(out, formula_calls(f))
    return(out)
  }
  if (inherits(x, c("frmtmb_formula", "brmsformula"))) {
    out <- add(out, response_formula_calls(x[["formula"]]))
    for (p in c(x[["pforms"]], x[["nlforms"]])) {
      out <- add(out, list(aterm = character(0),
                           other = expr_call_names(p)))
    }
    return(out)
  }
  # a bare formula reaching here is a response formula
  response_formula_calls(x)
}

#' Split ONE response formula at the bar on its left-hand side.
#'
#' @noRd
response_formula_calls <- function(f) {
  if (!inherits(f, "formula") || length(f) < 3L) {
    return(list(aterm = character(0), other = expr_call_names(f)))
  }
  lhs <- f[[2L]]
  rest <- expr_call_names(f[[3L]])
  if (is.call(lhs) && is.name(lhs[[1L]]) &&
        identical(as.character(lhs[[1L]]), "|")) {
    # y | <aterms> ~ rhs. The response expression itself is lhs[[2]] and
    # is not aterm position; cbind() and friends live there.
    list(aterm = expr_call_names(lhs[[3L]]),
         other = unique(c(expr_call_names(lhs[[2L]]), rest)))
  } else {
    list(aterm = character(0),
         other = unique(c(expr_call_names(lhs), rest)))
  }
}

#' Every function called anywhere in one expression.
#'
#' `pkg::fn` resolves to `fn`, which is deliberate: a qualified call is
#' the same feature. It is also why a THIRD party's `somepkg::s(x)`
#' matches `s()`, which is one of the holes the header comment names.
#'
#' @noRd
expr_call_names <- function(x, out = character(0)) {
  if (inherits(x, "formula") || is.call(x)) {
    if (!inherits(x, "formula")) {
      fn <- x[[1L]]
      nm <- if (is.name(fn)) {
        as.character(fn)
      } else if (is.call(fn) &&
                   identical(as.character(fn[[1L]]), "::")) {
        as.character(fn[[3L]])
      } else {
        ""
      }
      if (nzchar(nm)) out <- c(out, nm)
    }
    for (i in seq_along(x)) {
      e <- x[[i]]
      # match.call() leaves an unsupplied argument as the empty symbol,
      # which is neither NULL nor missing to is.null()
      if (is.name(e) && !nzchar(as.character(e))) next
      out <- expr_call_names(e, out)
    }
  }
  unique(out)
}

#' The family names of a model, from whichever of the two routes it
#' arrived on.
#'
#' It does NOT decompose a mixture. `mixture(gaussian, gaussian)` has
#' the single family name `"mixture(gaussian, gaussian)"`, on the
#' `bf()` object and on the fitted object alike, so a refused
#' `gaussian` row misses a mixture of gaussians while reaching a plain
#' gaussian correctly. A miss, not a false alarm.
#'
#' @noRd
model_family_names <- function(x, family = NULL) {
  nm <- function(f) {
    if (is.null(f)) return(character(0))
    if (is.character(f)) return(f)
    if (is.function(f)) f <- tryCatch(f(), error = function(e) NULL)
    as.character(f[["family"]] %||% character(0))
  }
  if (inherits(x, "frmtmb_fit")) {
    return(unique(unlist(lapply(x[["spec"]][["responses"]],
                                function(r) nm(r[["family"]])))))
  }
  bf_fams <- if (inherits(x, "frmtmb_mvformula")) {
    unlist(lapply(x[["forms"]], function(f) nm(f[["family"]])))
  } else if (inherits(x, "frmtmb_formula")) {
    nm(x[["family"]])
  }
  unique(c(nm(family), bf_fams))
}

#' Is this name a function of a base-priority R package?
#'
#' `trunc()`, `weights()` and `ar()` are registry display names AND
#' functions of base or stats, so a formula that writes one of them is
#' as likely to be doing arithmetic as to be using the feature. All
#' fourteen base-priority packages are checked, because they are the
#' ones present in every R installation; that is what makes the answer
#' independent of what the session has attached, which a search-path
#' test would not be. `inherits = FALSE`, so a user's own masking
#' definition cannot flip it either way.
#'
#' Memoized, because a name that matches nothing reads all fourteen
#' packages to find that out. The answer cannot change within a
#' session: base-priority packages ship with R and are not replaced.
#'
#' A package whose namespace is not loaded is not loaded to answer:
#' this runs on a path deliberately cleared of package loading, and
#' loading tcltk on a headless Linux machine warns "no DISPLAY variable
#' so Tk is not available" into a call that has nothing to do with Tk.
#' That warning escaped the Ubuntu check of frmtmb.sample on d24f7b86.
#' Such a package is read from its lazy-load database
#' (`base_r_objects()`) and decided by `exists()` alone: reading the
#' value would unserialize a closure whose environment is the package
#' namespace, and R resolves that reference by loading the package
#' (`tclVar` loads tcltk, `mle` loads stats4). The cost is that a
#' non-function object of such a package counts as a function; none of
#' the registry's call-shaped names is one (`dev/ciharden-regexists.R`).
#'
#' @noRd
base_r_function <- local({
  seen <- new.env(parent = emptyenv())
  function(nm) {
    if (!is.null(seen[[nm]])) return(seen[[nm]])
    pkgs <- c("base", "stats", "utils", "methods", "graphics",
              "grDevices", "splines", "stats4", "grid", "tools",
              "parallel", "compiler", "datasets", "tcltk")
    ans <- FALSE
    for (p in pkgs) {
      ns <- base_r_objects(p)
      if (!is.null(ns) && exists(nm, ns, inherits = FALSE) &&
            (!isNamespace(ns) || is.function(get(nm, ns)))) {
        ans <- TRUE
        break
      }
    }
    seen[[nm]] <- ans
    ans
  }
})

#' The objects of a base-priority package without loading it: its
#' namespace when that is loaded, otherwise its lazy-load database read
#' into an environment of promises, which runs no `.onLoad()` as long as
#' no promise is forced (`base_r_function()` only asks `exists()`). NULL
#' when the package has no database, as `datasets` has none of code.
#'
#' @noRd
base_r_objects <- local({
  dbs <- new.env(parent = emptyenv())
  function(p) {
    if (isNamespaceLoaded(p)) return(asNamespace(p))
    if (!is.null(dbs[[p]])) return(dbs[[p]])
    dir <- find.package(p, lib.loc = .Library, quiet = TRUE)
    base <- file.path(dir, "R", p)
    if (!length(dir) || !file.exists(paste0(base, ".rdx"))) return(NULL)
    e <- new.env(parent = emptyenv())
    ok <- tryCatch({
      lazyLoad(base, envir = e)
      TRUE
    }, error = function(err) FALSE)
    if (!ok) return(NULL)
    dbs[[p]] <- e
    e
  }
})

#' The name a formula writes for a call-shaped feature.
#'
#' The DISPLAY name minus its parentheses, not the registry key. The
#' key is the parser name and reads better, but three keys are shared
#' between two features each, so keying on it forced a downgrade to a
#' warning on a shared key, and that downgrade could be triggered by
#' anybody's later registration. See the header comment.
#'
#' @noRd
feature_call_name <- function(feat) sub("[(][)]$", "", feat)

#' What a refused row can be matched against, and why not, when not.
#'
#' Returns `"family"`, `"aterm"`, `"call"`, or a sentence naming the
#' reason the row cannot be decided from a formula. Everything it reads
#' comes from the registry and from the installed base packages at call
#' time, so a collision registered later is detected rather than
#' assumed away.
#'
#' @noRd
preflight_route <- function(feat, kind, ft) {
  if (identical(kind, "family")) return("family")
  if (!grepl("[(][)]$", feat)) {
    return(paste0("the vocabulary writes it without parentheses, so ",
                  "it is not declared as a call a formula can hold"))
  }
  # An addition term has a position of its own and needs no name test.
  if (identical(kind, "aterm")) return("aterm")
  nm <- feature_call_name(feat)
  # No uniqueness test is needed, and none is written, because two
  # call-shaped rows cannot want the same call name. Display names are
  # unique (core refuses a second registration of one: "One display
  # name carries one kind"), dropping a fixed "()" suffix is injective
  # over the names that carry it, and a bare name cannot reach this
  # point at all because the test above returns first. That is what
  # makes a registered refusal un-revocable; see the header comment.
  #
  # The display name is not always what a formula writes, though, and
  # where it is not, matching it would MISS IN SILENCE while every
  # other undecidable case warns. The registry says which those are:
  # its key is the parser name, so key and display name agree for every
  # call-shaped row except `mi_pred()`, `gp_pred()` and `cs_pred()`,
  # whose keys are `mi`, `gp` and `cs`. Those three warn.
  key <- ft$key[match(feat, ft$name)]
  if (!is.na(key) && !identical(key, nm)) {
    return(paste0("no formula writes the call '", nm, "': the ",
                  "vocabulary spells this feature '", feat,
                  "' and a formula writes it '", key,
                  "', which another feature also writes"))
  }
  if (base_r_function(nm)) {
    return(paste0("'", nm, "' is also a base R function, so a call to ",
                  "it in a formula need not be the feature at all"))
  }
  "call"
}

#' Stop when the registry declares this sampler refused for something
#' the model uses.
#'
#' Runs before anything is taped, before the formula route calls
#' `frm()`, and before the rstan and tmbstan namespaces are loaded, so
#' the cost of a refused model is the registry read and nothing else.
#'
#' `rows` defaults to the registry and is an argument so that the tests
#' can drive each branch of the decision. Registering a synthetic
#' feature instead would append to a registry that lives as long as the
#' session, and every later test file would inherit it.
#'
#' @noRd
sample_preflight <- function(x, family = NULL, rows = NULL,
                             what = "frm_sample()") {
  rows <- rows %||% frmtmb::frm_compat("frm_sample", status = "refused")
  if (!nrow(rows)) return(invisible(NULL))
  ft <- frmtmb::frm_compat_features()
  calls <- formula_calls(
    if (inherits(x, "frmtmb_fit")) x[["bform"]] else x)
  fams <- model_family_names(x, family)
  for (i in seq_len(nrow(rows))) {
    feat <- rows$feature_b[[i]]
    route <- preflight_route(feat, rows$kind_b[[i]], ft)
    nm <- feature_call_name(feat)
    hit <- switch(route,
                  family = nm %in% fams,
                  aterm = nm %in% calls$aterm,
                  call = nm %in% calls$other,
                  {
                    # Undecidable from a formula. Saying nothing would
                    # make a registered refusal look like a passed
                    # check; refusing anyway would fire on models that
                    # do not use the feature at all.
                    frm_warning(what, ": the registry refuses this sampler ",
                                "with '", feat, "', and the pre-flight ",
                                "cannot tell whether this model uses it: ",
                                route, ". Sampling anyway. See ",
                                "frm_compat(\"frm_sample\", \"", feat,
                                "\")", call. = FALSE)
                    FALSE
                  })
    if (isTRUE(hit)) {
      frm_stop(what, " is refused on a model using ", feat, ": ",
               rows$note[[i]], " (registry row: ",
               "frm_compat(\"frm_sample\", \"", feat, "\"))",
               call. = FALSE)
    }
  }
  invisible(NULL)
}

# --- sampler start values and bounds -----------------------------------

#' Pull a start value strictly inside the bounding box. Stan turns a
#' bound into a constrained transform, and a start AT the bound maps to
#' an infinite unconstrained value: rstan then reports "Initialization
#' failed" and names neither the parameter nor the bound. The pad is
#' relative so it survives bounds of any magnitude; a box narrower than
#' two pads collapses to its midpoint.
#'
#' @noRd
clamp_into_bounds <- function(v, lower, upper) {
  if (is.null(lower) && is.null(upper)) return(v)
  lo <- lower %||% rep(-Inf, length(v))
  hi <- upper %||% rep(Inf, length(v))
  pad_lo <- 1e-3 * pmax(1, abs(lo))
  pad_hi <- 1e-3 * pmax(1, abs(hi))
  fl <- ifelse(is.finite(lo), lo + pad_lo, -Inf)
  ce <- ifelse(is.finite(hi), hi - pad_hi, Inf)
  narrow <- is.finite(lo) & is.finite(hi) & fl > ce
  mid <- (lo + hi) / 2
  fl[narrow] <- mid[narrow]
  ce[narrow] <- mid[narrow]
  pmin(pmax(v, fl), ce)
}

#' Per-chain initial values around the ML mode: chain 1 exactly at the
#' mode (the short-warmup anchor), later chains at mode + N(0, jitter)
#' on the unconstrained scale, restoring the overdispersion Rhat needs
#' to detect chains agreeing for the wrong reason. Every chain - the
#' mode-anchored one included, because the bounds can exclude the ML
#' mode itself - is then pulled inside the bounds.
#'
#' @noRd
mode_inits <- function(mode, chains, jitter, lower = NULL, upper = NULL) {
  mode <- as.numeric(mode)
  if (!is.finite(jitter) || jitter <= 0 || chains <= 1L) {
    v <- clamp_into_bounds(mode, lower, upper)
    return(lapply(seq_len(max(chains, 1L)), function(i) v))
  }
  lapply(seq_len(chains), function(i) {
    v <- if (i == 1L) mode else mode + stats::rnorm(length(mode), 0, jitter)
    clamp_into_bounds(v, lower, upper)
  })
}

#' Evaluate `expr` on R's random stream seeded by `seed`, and leave the
#' caller's stream as it was. brms's `seed` reaches Stan alone, because
#' Stan draws brms's inits; frm_sample()'s jitter is drawn in R, and a
#' seed that fixed Stan and not the jitter left two seeded calls with
#' different chains. `NULL` or `NA`, brms's "no seed", draws from the
#' session's stream, so `set.seed()` still decides.
#'
#' @noRd
with_sample_seed <- function(seed, expr) {
  if (is.null(seed) || (length(seed) == 1L && is.na(seed))) return(expr)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    frm_stop("frm_sample(seed = ) takes one finite number, or NA for no ",
             "seed; got ", deparse1(seed), call. = FALSE)
  }
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  })
  set.seed(seed)
  expr
}

#' Per-chain numeric inits carried as one-dimensional arrays.
#'
#' tmbstan's model has one parameter, `vector[N] y`, and hands a
#' numeric init to rstan as `y`. rstan reads a length-one numeric as a
#' SCALAR, so on a model with one outer parameter (`y ~ 0 + x`) every
#' chain died with "no more scalars to read" before sampling. A
#' `dim` attribute makes it the vector of length one Stan declares, and
#' tmbstan's own init sanitizer keeps the attribute.
#'
#' @noRd
stan_init_arrays <- function(init) {
  if (!is.list(init)) return(init)
  lapply(init, function(v) {
    if (is.numeric(v) && is.null(dim(v))) array(v, dim = length(v)) else v
  })
}

#' tmbstan widens outer-length bounds over the whole parameter vector
#' when the objective has random effects (the inner block is unbounded);
#' the inits are the full vector too, so they must be clamped against
#' the same widening.
#'
#' @noRd
mode_aligned_bounds <- function(obj, bounds, laplace, n) {
  if (is.null(bounds)) return(NULL)
  if (length(bounds$lower) == n) return(bounds)
  rnd <- obj$env$random
  if (laplace || is.null(rnd) || length(rnd) + length(bounds$lower) != n) {
    return(NULL)
  }
  lo <- rep(-Inf, n)
  hi <- rep(Inf, n)
  lo[-rnd] <- bounds$lower
  hi[-rnd] <- bounds$upper
  list(lower = lo, upper = hi)
}

#' Retape the fit's objective with priors added; parameters start at the
#' ML estimates so sampling initializes at (near) the posterior mode.
#'
#' A model with `mo()` terms is retaped in softmax coordinates with
#' brms's `dirichlet(1)` on each simplex (`mo_sample_frame()`); every
#' other model is retaped as it always was.
#'
#' @noRd
prior_augmented_obj <- function(fit, entries, chart = NULL) {
  nll <- build_objective(mo_sample_frame(fit))
  entries <- mo_sample_entries(entries, fit)
  nlp <- neg_log_prior_fn(entries)
  snlp <- sample_extra_nlp(fit, entries)
  tpl <- fit$frame[["par_template"]]
  # [[ ]] to avoid $ partial matching ("b" matching "beta" in GLMs)
  random <- c(if (!is.null(tpl[["b"]])) "b",
              if (!is.null(tpl[["miss"]])) "miss")
  if (!length(random)) random <- NULL
  fn <- if (is.null(snlp)) {
    function(pars) nll(pars) + nlp(pars)
  } else {
    function(pars) nll(pars) + nlp(pars) + snlp(pars)
  }
  RTMB::MakeADFun(chart_fn(fn, chart),
                  chart_pars(mo_sample_pars(fit, fit$estimates), chart),
                  random = random, map = fit$frame[["map"]], silent = TRUE)
}

#' The sampling objective read in the charts of core's
#' `autocor_bound_chart()`: each bounded residual-correlation parameter
#' is sampled as `v`, one logistic onto its natural interval, and handed
#' to the model as the internal parameter, with `log |d internal / d v|`
#' added to the density. Stan's log bound on the internal parameter,
#' stacked on the term's own map, gave 5 or 6 divergent transitions per
#' 4 chains of 2000 on a cosy near 0, where brms and this chart have
#' none (dev/acprior-cosydiv.R).
#'
#' @noRd
chart_fn <- function(fn, chart) {
  if (!length(chart)) return(fn)
  function(pars) {
    "[<-" <- RTMB::ADoverload("[<-")
    jac <- 0
    for (tm in chart) {
      v <- pars[[tm$comp]][tm$idx]
      x <- pars[[tm$comp]]
      x[tm$idx] <- tm$to_internal(v)
      pars[[tm$comp]] <- x
      jac <- jac + tm$logjac(v)
    }
    fn(pars) - jac
  }
}

#' A parameter list moved into those charts, for the start.
#'
#' @noRd
chart_pars <- function(pars, chart) {
  for (tm in chart) {
    pars[[tm$comp]][tm$idx] <- tm$from_internal(pars[[tm$comp]][tm$idx])
  }
  pars
}

#' A warning for each charted parameter whose ML estimate lies outside
#' its sampled interval, so the chains cannot start at the mode. When
#' the bound is the default's (cosy's 0, brms's declaration) the user
#' asked for nothing, so the message says where the bound comes from and
#' how to sample frm()'s whole window instead.
#'
#' @noRd
chart_mode_warning <- function(fit, chart) {
  if (!length(chart)) return(invisible(NULL))
  ap <- autocor_flat_priors(fit)
  multi <- length(fit$spec$responses) > 1L
  th <- fit$estimates[["thetaac"]]
  for (tm in chart) {
    est <- tm$to_natural(th[tm$idx])
    if (est > tm$L && est < tm$U) next
    dflt <- ap$class == tm$class & ap$resp == tm$resp & !is.na(ap$lb) &
      ap$lb == tm$L
    rs <- if (multi) paste0(", resp = \"", tm$resp, "\"") else ""
    if (any(dflt) && est <= tm$L) {
      frm_warning("The ML estimate of ", tm$class, ", ",
                  format(est, digits = 4), ", is below ", tm$L, ", the ",
                  "lower bound brms declares for ", tm$class, ", which is ",
                  "frm_sample()'s default there; every chain starts just ",
                  "inside it. To sample the whole window frm() fits, write ",
                  "set_prior(\"\", class = \"", tm$class, "\"", rs,
                  ", lb = -1/", round(-1 / tm$image[1L]), ")",
                  call. = FALSE)
    } else {
      frm_warning("The ML estimate of ", tm$class, ", ",
                  format(est, digits = 4), ", is outside the requested ",
                  "interval (", format(tm$L, digits = 6), ", ",
                  format(tm$U, digits = 6), "); every chain starts just ",
                  "inside it. The bounded posterior is not centered on the ",
                  "unconstrained ML estimate", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Draws moved back from the charts to the internal parameters, which
#' everything after the sampler reads.
#'
#' @noRd
chart_draws_back <- function(m, fit, chart, laplace) {
  if (!length(chart)) return(m)
  pidx <- draws_par_index(fit, laplace = laplace)
  for (tm in chart) {
    j <- pidx[[tm$comp]][tm$idx]
    if (is.null(j) || is.na(j) || j > ncol(m)) next
    m[, j] <- tm$to_internal(m[, j])
  }
  m
}

#' The densities the sampler adds outside the resolved entries: brms's
#' `dirichlet(1)` on each `mo()` simplex and on a mixture's weights.
#'
#' @noRd
sample_extra_nlp <- function(fit, entries) {
  f1 <- mo_simplex_nlp(fit, entries)
  f2 <- mix_theta_nlp(fit)
  if (is.null(f1)) return(f2)
  if (is.null(f2)) return(f1)
  function(pars) f1(pars) + f2(pars)
}

# --- mo() simplexes ---------------------------------------------------
#
# The fit holds a mo() simplex in the chart of frmtmb's mo_simplex(),
# whose faces are folds at finite coordinates: right for maximum
# likelihood, wrong for a sampler. With no density on the coordinates
# the chart itself is the prior, prod(w)^(-1/2) summed over 2^D sign
# sheets, which is not integrable at the barycenter: a weakly
# identified simplex sat there with sd 0 and coordinates near 1e18
# (dev/optima-rev-sample2-arm.R). brms's dirichlet(1) carried to that
# chart vanishes on every fold, so its 2^D sheets are walled off from
# each other: 50 divergences, chains stuck on different sheets (R-hat
# 8.47). The softmax is one sheet with no fold, and the same prior
# there samples cleanly. So the sampler runs in softmax coordinates
# and the draws are carried back to the fit's chart before anything
# reads them (dev/optima-findings.md, punch round 1).

#' The frame the sampler tapes: the fit's own, or a copy whose `mo()`
#' simplexes read softmax coordinates.
#'
#' @noRd
mo_sample_frame <- function(fit) {
  fr <- if (!length(simplex_frame_terms(fit))) fit$frame else {
    mo_chart_frame(fit$frame, "softmax")
  }
  # a scalar GEV xi in brms's tmp_xi chart (frmtmb's gev_xi_chart_frame())
  gev_xi_chart_frame(fr, gev_xi_chart_terms(fit))
}

#' A parameter list with each `mo()` simplex moved from the fit's chart
#' to softmax coordinates (a weight at a face starts at 1e-8, since a
#' face is at an infinite softmax coordinate).
#'
#' @noRd
mo_sample_pars <- function(fit, pars) {
  for (tm in simplex_frame_terms(fit)) {
    pars[[tm$zeta]] <- mo_coords(mo_simplex(pars[[tm$zeta]]), "softmax")
  }
  gev_xi_chart_start(fit, gev_xi_chart_terms(fit), pars)
}

#' The `mo()` simplexes (their `zeta` components) a fit's own prior puts
#' a density on: a MAP fit with a `simo` row. Empty for an ML fit.
#'
#' @noRd
mo_fit_priored <- function(fit) {
  if (is.null(fit$prior) || !length(simplex_frame_terms(fit))) {
    return(character(0))
  }
  ent <- tryCatch(resolve_prior_input(fit, fit$prior)$entries,
                  error = function(e) list())
  unlist(lapply(ent, function(e) {
    if (identical(e$dist$kind, "dirichlet")) e$comp
  }))
}

#' Resolved prior entries for the sampler's tape: a Dirichlet on a `mo()`
#' simplex (class `"simo"`) is read in the softmax chart the sampler
#' retapes the simplex in (`mo_sample_frame()`), with that chart's
#' Jacobian (frmtmb's `dirichlet_logdens()`), where the fit's entry was
#' resolved in the fit's own chart.
#'
#' @noRd
mo_sample_entries <- function(entries, fit = NULL) {
  # without the fit there is no GEV chart to read a prior in
  terms <- if (!is.null(fit)) {
    Filter(function(tm) !length(tm$blocked), gev_xi_chart_terms(fit))
  }
  lapply(entries, function(e) {
    if (identical(e$dist$kind, "dirichlet")) e$dist$chart <- "softmax"
    # a prior on a scalar xi is brms's prior on tmp_xi, the coefficient
    # the chart samples, so it is read on that coefficient as it stands
    for (tm in terms) {
      if (identical(e$comp, tm$par) && length(e$idx) &&
            all(e$idx %in% tm$idx)) {
        e$scale <- "tmp_xi"
        e$link <- NULL
      }
    }
    e
  })
}

#' brms's default prior on a `mo()` simplex, `dirichlet(1)`, the uniform
#' density on the simplex, as a negative log density of the softmax
#' coordinates: the uniform density times the softmax's Jacobian,
#' `prod(w)` over all D weights, is `Gamma(D) prod(w)`. `NULL` without
#' `mo()` terms.
#'
#' It is not a prior frm_sample(prior = "flat") turns off: a flat
#' density on the simplex IS dirichlet(1), and a flat density on the
#' coordinates instead is `1 / prod(w)` on the simplex, infinite mass at
#' every face (the softmax arm of dev/optima-rev-sample2-arm.R: every
#' draw at a vertex, 44 divergences).
#'
#' A simplex with a `simo` prior of its own (a resolved Dirichlet entry,
#' `mo_sample_entries()`) takes that density instead, so it is left out
#' here rather than given two.
#'
#' @noRd
mo_simplex_nlp <- function(fit, entries = list()) {
  priored <- unlist(lapply(entries, function(e) {
    if (identical(e$dist$kind, "dirichlet")) e$comp
  }))
  terms <- Filter(function(tm) !tm$zeta %in% priored,
                  simplex_frame_terms(fit))
  if (!length(terms)) return(NULL)
  function(pars) {
    "c" <- RTMB::ADoverload("c")
    out <- 0
    for (tm in terms) {
      a <- c(0, pars[[tm$zeta]])
      # log(sum(exp(a))) without overflow on the tape
      lse <- a[1L]
      for (k in seq_along(a)[-1L]) lse <- RTMB::logspace_add(lse, a[k])
      # sum(log w) = sum(a) - D * lse
      out <- out - (sum(a) - tm$D * lse) - lgamma(tm$D)
    }
    out
  }
}

#' Why a scalar GEV `xi` is not sampled through `tmp_xi`, and what that
#' costs, for the warning and the disclosure note.
#'
#' @noRd
gev_xi_blocked_text <- function(tm, of = "") {
  paste0(tm$dpar %||% "xi", of, " of gen_extreme_value() is sampled ",
         "on its ", tm$link$name, " scale, not through brms's tmp_xi, ",
         "because of ", paste(tm$blocked, collapse = " and "), ", and ",
         "a prior on it is a density on xi itself. The support's edge is ",
         "then a wall the sampler meets: on the data of ",
         "dev/famuni-out/tmpxi.txt, 957 to 3388 of 4000 transitions ",
         "diverged that way and none through tmp_xi. Check the divergences")
}

#' The draws matrix with a scalar GEV `xi` carried from `tmp_xi`, the
#' chart the sampler ran in, to the fit's own coefficient: each draw's
#' `xi` is the one its density used, which the tape REPORTs, under the
#' fit's link. brms's `compute_xi()` does the same after sampling.
#'
#' @noRd
gev_draws_to_xi <- function(m, obj, fit, terms, laplace) {
  terms <- Filter(function(tm) !length(tm$blocked), terms)
  if (!length(terms)) return(m)
  tm <- terms[[1L]]
  pidx <- draws_par_index(fit, laplace = laplace)
  tpl <- fit$frame[["par_template"]]
  fx <- if (identical(tm$par, "betad")) fit$frame[["betad_fixed_idx"]]
  free <- setdiff(seq_along(tpl[[tm$par]]), fx)
  j <- pidx[[tm$par]][match(tm$idx, free)]
  np <- if (laplace) length(obj$par) else length(obj$env$last.par)
  for (i in seq_len(nrow(m))) {
    row <- as.numeric(m[i, seq_len(np)])
    rep <- if (laplace) {
      obj$fn(row)
      obj$report()
    } else {
      obj$report(row)
    }
    m[i, j] <- tm$link$linkfun(rep[["frmtmb_gev_xi"]][1L])
  }
  m
}

#' The draws matrix with each `mo()` simplex's softmax columns carried to
#' the fit's chart, the coordinates every reader of the draws hands back
#' to the model (`draws_to_natural()` names them as brms's `simo_`).
#'
#' @noRd
mo_draws_to_chart <- function(m, fit, laplace) {
  terms <- simplex_frame_terms(fit)
  if (!length(terms)) return(m)
  pidx <- draws_par_index(fit, laplace = laplace)
  for (tm in terms) {
    j <- pidx[[tm$zeta]]
    if (is.null(j) || max(j) > ncol(m)) next
    m[, j] <- draws_rowmap(m[, j, drop = FALSE], function(z) {
      mo_coords(mo_simplex(z, "softmax"))
    })
  }
  m
}

#' `draws_to_natural()`'s blocks for the `mo()` simplexes, in the shape
#' of an ordinal block (`draws_ordinal_block()`): the `D - 1` chart
#' columns become brms's `simo_<label>1[1..D-1]` in place, and
#' `simo_<label>1[D]` is added right after them (`adjacent`), so each
#' term's weights sit together in brms's order.
#'
#' @noRd
draws_mo_cols <- function(fit) {
  lapply(simplex_frame_terms(fit), function(tm) {
    list(internal = paste0(tm$zeta, "_", seq_len(tm$D - 1L)),
         names = tm$names,
         map = function(r) mo_simplex(r),
         inv = function(v) mo_coords(v[seq_len(tm$D)]),
         adjacent = TRUE)
  })
}

# --- non-centered sampling --------------------------------------------
#
# See R/covstruct.R for the factor accessors and for which structures
# have one. Everything here is SAMPLING-ONLY: the fit, its tape and every
# number it reports are untouched, because the Laplace approximation
# integrates `b` out and that integral does not care how the integrated
# variable is parameterized.

#' Why a structure has no non-centered form. Read out loud by the
#' message that names the blocks left centered, so a fallback is never
#' silent.
#'
#' @noRd
ncp_reasons <- c(
  us_t = "Student-t latent: a scale mixture, not a linear factor",
  diag_t = "Student-t latent: a scale mixture, not a linear factor",
  car = "sparse CAR precision; its factor is dense",
  spde = "sparse SPDE precision; its factor is dense",
  gr_prec = "sparse precision from gr(prec = ); its factor is dense",
  ou = "dense covariance over the field's positions",
  exp = "dense covariance over the field's positions",
  gau = "dense covariance over the field's positions",
  mat = "dense covariance over the field's positions",
  toep = "not positive definite everywhere in its parameterization",
  homtoep = "not positive definite everywhere in its parameterization"
)

#' Why this block stays centered.
#'
#' @noRd
ncp_reason <- function(bk) {
  cs <- bk[["covstruct"]]
  if (is_student_block(bk)) return(unname(ncp_reasons[["us_t"]]))
  # [[ ]] on a named character vector ERRORS on a missing name, so the
  # membership test comes first
  if (cs %in% names(ncp_reasons)) return(unname(ncp_reasons[[cs]]))
  if (!covstruct_has_chol(bk)) {
    return("no Cholesky factor is registered for this structure")
  }
  # a factor exists, but some parameter of the block is neither a
  # standard deviation nor a correlation the prior lane can name, so
  # nothing would bound the direction the transform opens
  paste("it has a parameter that is neither a standard deviation nor a",
        "correlation with a density on it, so a non-centered chain",
        "would be free in a direction no prior bounds")
}

#' The `theta` positions a resolved prior specification covers.
#'
#' @noRd
ncp_priored_theta <- function(entries) {
  unlist(lapply(entries %||% list(), function(e) {
    if (identical(e$comp, "theta")) e$idx
  })) %||% integer(0)
}

#' Which blocks a sampling run non-centers, and which it leaves alone.
#'
#' TWO conditions, and they are the same condition twice. A block is
#' non-centered only when (a) every parameter it has is a standard
#' deviation or a correlation, with a registered factor (see
#' R/covstruct.R), and (b) every one of those parameters carries a
#' PRIOR. Both are about not handing the sampler a direction it can run
#' away in.
#'
#' Non-centering gives NUTS the run of theta's whole range. Under a FLAT
#' prior on a log standard deviation that range has a flat tail at the
#' bottom: send `sd` down and `z` stays where it is, `b = sd z` goes to
#' zero, the likelihood settles on the no-random-effect model and the
#' density stops changing. The centered funnel is the only thing that
#' keeps a chain out of it. Measured on a six-group random-intercept
#' model, one chain of 2000, three seeds: sampled from the FIT, where
#' every prior is flat, the non-centered chain walks `theta` to -1e15
#' with a bulk-ESS of 1; sampled from the FORMULA, where the default
#' `student_t(3, 0, s)` on the sd makes that tail integrable, it runs at
#' a bulk-ESS of 174-284 against 3-48 centered. Adding the prior by hand
#' on the fit route does the same thing there (297-355).
#'
#' A CORRELATED block is the same argument with a different parameter.
#' Flat on frmtmb's correlation parameter is `(1 - rho^2)^-3/2`, which
#' is improper, and before 0.39 that kept every correlated block
#' centered whatever else was priored. The LKJ prior (R/priors.R) closes
#' that tail, the formula route applies `lkj(1)` by default, and the
#' rule above then admits those blocks with nothing else changed.
#'
#' `laplace = TRUE` integrates the random effects out, so there is no
#' funnel left to fix and the transform would only cost inner Newton
#' steps: that route stays centered whatever was asked for.
#'
#' @noRd
ncp_plan <- function(fit, reparameterize, laplace, entries = NULL) {
  blocks <- fit$frame[["re_blocks"]] %||% list()
  if (!isTRUE(reparameterize) || !length(blocks) || laplace) {
    return(list(idx = integer(0), centered = character(0),
                labels = character(0)))
  }
  priored <- ncp_priored_theta(entries)
  has_prior <- vapply(blocks, function(bk) {
    all(bk[["theta_idx"]] %in% priored)
  }, TRUE)
  ok <- vapply(blocks, ncp_eligible, TRUE) & has_prior
  # an rr() block's own coefficients ARE standard normal factors, so it
  # is non-centered already and has nothing to report
  done <- vapply(blocks, function(bk) {
    identical(bk[["covstruct"]], "rr")
  }, TRUE)
  lab <- vapply(blocks, function(bk) {
    paste0(bk[["term_label"]], " [", bk[["covstruct"]], "]")
  }, "")
  left <- which(!ok & !done)
  list(idx = which(ok), labels = unname(lab[ok]),
       centered = unname(vapply(left, function(i) {
         bk <- blocks[[i]]
         why <- if (!ncp_eligible(bk)) {
           ncp_reason(bk)
         } else if (block_n_cor(bk) > 0L) {
           # the correlation is the interesting half here: flat on
           # frmtmb's correlation parameter is (1 - rho^2)^-3/2, which
           # is improper, and the centered geometry was what kept a
           # chain out of that tail
           paste("its variance or correlation parameters have a flat",
                 "prior here, and flat on a correlation is",
                 "(1 - rho^2)^-3/2, improper. Give the block priors,",
                 "set_prior(class = \"sd\") and",
                 "set_prior(class = \"cor\"), which frm_sample()",
                 "supplies by default unless prior = \"flat\" turned",
                 "them off")
         } else if (identical(bk[["covstruct"]], "hsgp")) {
           # the sd advice below would leave the lengthscales flat
           paste("its sd or its length scales have a flat prior here.",
                 "Give them priors, set_prior(class = \"sdgp\") and",
                 "set_prior(class = \"lscale\"), which frm_sample()",
                 "supplies by default unless prior = \"flat\" turned",
                 "them off")
         } else {
           paste("its variance parameter has a flat prior here, and a",
                 "non-centered chain walks the flat tail that opens at",
                 "sd = 0. Give it a prior, set_prior(class = \"sd\"),",
                 "which frm_sample() supplies by default unless",
                 "prior = \"flat\" turned it off")
         }
         paste0(lab[i], ": ", why)
       }, "")))
}

#' One line per block that could not be non-centered, so a partial
#' fallback is visible without being an error. `message()`, not
#' `warning()`: the run is correct either way, only slower to mix.
#'
#' @noRd
announce_ncp <- function(plan) {
  if (!length(plan$centered)) return(invisible(NULL))
  hdr <- if (!length(plan$idx)) {
    paste("frm_sample(): sampling stays centered: no random-effect",
          "block of this model has a non-centered form:")
  } else {
    paste("frm_sample(): non-centered where possible; these blocks",
          "stay centered:")
  }
  frm_message(paste(c(hdr, paste0("  ", plan$centered)), collapse = "\n"))
  invisible(NULL)
}

#' The sampling objective with the planned blocks reparameterized, and
#' the prior terms (if any) added the way prior_augmented_obj() adds
#' them. The frame copy carrying `ncp_blocks` never leaves this
#' function, so no fitted object can acquire a non-centered tape.
#'
#' @noRd
ncp_objective <- function(fit, idx, entries, chart = NULL) {
  frame <- mo_sample_frame(fit)
  frame[["ncp_blocks"]] <- idx
  nll <- build_objective(frame)
  entries <- mo_sample_entries(entries, fit)
  fn <- if (length(entries)) {
    nlp <- neg_log_prior_fn(entries)
    function(pars) nll(pars) + nlp(pars)
  } else {
    nll
  }
  snlp <- sample_extra_nlp(fit, entries)
  if (!is.null(snlp)) {
    fn0 <- fn
    fn <- function(pars) fn0(pars) + snlp(pars)
  }
  tpl <- fit$frame[["par_template"]]
  # [[ ]] to avoid $ partial matching ("b" matching "beta" in GLMs)
  random <- c(if (!is.null(tpl[["b"]])) "b",
              if (!is.null(tpl[["miss"]])) "miss")
  if (!length(random)) random <- NULL
  RTMB::MakeADFun(chart_fn(fn, chart),
                  chart_pars(mo_sample_pars(fit, ncp_start_pars(fit, idx)),
                             chart),
                  random = random, map = fit$frame[["map"]], silent = TRUE)
}

#' The starting parameter list on the non-centered scale: `z0 = L^-1 b`
#' per reparameterized block, so that `init = "last.par.best"` still
#' starts the chains at the ML mode of the MODEL rather than at the
#' mode of one parameterization and the origin of the other.
#'
#' @noRd
ncp_start_pars <- function(fit, idx) {
  pars <- fit$estimates
  th <- pars[["theta"]] %||% numeric(0)
  bv <- pars[["b"]]
  for (i in idx) {
    bk <- fit$frame[["re_blocks"]][[i]]
    bv[bk[["b_idx"]]] <- ncp_unscale_b(bk,
                                       bv[bk[["b_idx"]]], th[bk[["theta_idx"]]])
  }
  pars[["b"]] <- bv
  pars
}

#' Put the draws matrix back on the centered scale, draw by draw: each
#' `z` maps through the `L(theta)` of ITS OWN draw. The random-effect
#' columns
#' that come out are the same quantity, in the same order and under the
#' same names, as the centered route's, which is what lets every draws
#' method downstream stay ignorant of the difference.
#'
#' @noRd
ncp_backtransform <- function(m, fit, idx) {
  pidx <- draws_par_index(fit)
  bc <- pidx[["b"]]
  tc <- pidx[["theta"]]
  blocks <- fit$frame[["re_blocks"]]
  for (r in seq_len(nrow(m))) {
    th <- unname(m[r, tc])
    for (i in idx) {
      bk <- blocks[[i]]
      cols <- bc[bk[["b_idx"]]]
      m[r, cols] <- ncp_scale_b(bk, unname(m[r, cols]),
                                th[bk[["theta_idx"]]])
    }
  }
  m
}

#' The columns of a draws matrix that are the model's outer parameters:
#' everything except the sampled random effects, the missing values and
#' `lp__`. `summary()`, `print()` and `check_laplace()` report these
#' alone, as brms's `summary()` leaves out `r_` and `lp__`.
#'
#' Decided by the fit's own labels rather than by a name pattern, so it
#' does not depend on which spelling a block's columns carry: the
#' `r_<group>[...]` of an ordinary block and the `b[<i>]` of a smooth or
#' reduced-rank one are both random.
#'
#' @noRd
draws_outer_cols <- function(x) {
  # a mixture's last weight is a derived column, and an outer parameter;
  # an ordinal block is stored under brms's names (draws_ordinal_cols())
  nc <- draws_natural_cols(x$fit)
  intersect(colnames(x$draws),
            c(brms_par_labels(x$fit, include_random = FALSE), nc$extra,
              unlist(lapply(nc$ordinal, `[[`, "names"))))
}

#' Whether a draws object came from `frm_sample(laplace = TRUE)`, which
#' samples the outer parameters alone. Counted rather than read off a
#' column name, because the random-effect columns carry brms's names
#' now and those depend on the model.
#'
#' @noRd
draws_is_laplace <- function(x) {
  fit <- x$fit
  n_full <- length(brms_par_labels(fit, include_random = TRUE))
  n_outer <- length(brms_par_labels(fit, include_random = FALSE))
  n_extra <- length(draws_natural_cols(fit)$extra)
  n_full > n_outer && ncol(x$draws) <= n_outer + 1L + n_extra
}

#' The core count tmbstan would actually use. rstan defaults its `cores`
#' argument to getOption("mc.cores", 1L), so a session-wide mc.cores
#' starts parallel chains just as an explicit cores = does, and the
#' Windows startup note has to read the same source or it never fires.
#'
#' @noRd
stan_cores <- function(args) {
  as.numeric(args$cores %||% getOption("mc.cores", 1L))[1L]
}

#' Packages this model's code comes from that a worker process could
#' not load. R serializes a closure's namespace by name, and when the
#' receiving process cannot load that namespace it substitutes the
#' global environment without stopping; the family then fails at its
#' first internal call with "could not find function". A namespace
#' built by pkgload::load_all() exists only in the session that built
#' it, so every parallel chain on Windows fails that way, and rstan
#' reports only that the chains returned no draws.
#'
#' @noRd
dev_namespaces_of <- function(fit) {
  fns <- list()
  for (resp in fit[["spec"]][["responses"]]) {
    fam <- resp[["family"]]
    fns <- c(fns, unlist(rapply(unclass(fam), function(f) f,
                                classes = "function", how = "list")))
  }
  nms <- vapply(fns, function(f) {
    e <- topenv(environment(f) %||% globalenv())
    if (isNamespace(e)) getNamespaceName(e) else ""
  }, "")
  nms <- unique(nms[nzchar(nms)])
  nms[vapply(nms, function(n) {
    exists(".__DEVTOOLS__", envir = asNamespace(n), inherits = FALSE)
  }, TRUE)]
}

#' Assemble the unfitted object the formula interface samples.
#'
#' @noRd
sample_assemble <- function(formula, data, family, data2, start,
                            fit_control, na.action, REML) {
  if (!inherits(formula, c("formula", "frmtmb_formula",
                           "frmtmb_mvformula"))) {
    frm_stop("frm_sample() takes a frmtmb fit or a formula; got an object ",
             "of class ", paste(class(formula), collapse = "/"),
             ". Fit with frm() first, or pass bf(y ~ x) with data =",
             call. = FALSE)
  }
  if (is.null(data)) {
    frm_stop("frm_sample() from a formula needs data =: there is no fitted ",
             "model to take the design from", call. = FALSE)
  }
  # no bounds here: the unfitted object is never optimized, and the box
  # Stan is given is resolved from the prior further down
  frm(formula, data, family = family, REML = REML, start = start,
      control = fit_control, na.action = na.action,
      data2 = data2, dry_run = "objective")
}

# --- default priors for frm_sample() -----------------------------------
#
# WHERE THEY APPLY. Both of frm_sample()'s routes, because both sample a
# POSTERIOR and the posterior is the whole answer on either. Until
# 0.43.0 the fit route was flat, on the argument that it exists to serve
# check_laplace(); the audit in dev/brms-vignette-audit.md measured what
# that costs a reader who ports a prior-free brm() call onto a fitted
# model, and it is not a diagnostic they get - it is the chain pinned at
# whatever boundary mode maximum likelihood found (kidney's sd(patient)
# at 3e-4), or a chain that mixes an order of magnitude worse than the
# same model sampled from the formula. check_laplace() keeps the
# unpenalized density it needs through frm_sample(.diagnostic = TRUE),
# which is now the only route that is flat by default.
#
# THE RECIPE is brms 2.23's, read off `brms::default_prior()` on matched
# models (tests/testthat/test-sample-direct.R checks it against brms
# directly where brms is installed). Writing y* for the response
# transformed by the mu link:
#
#   b (slopes)      flat, as brms leaves them
#   Intercept       student_t(3, round(median(y*), 1) - mean(offset), s),
#                   the offset being the predictor's own; brms takes it
#                   off only where y* was transformed at all
#   sd              student_t(3, 0, s) on the natural sd scale, written
#                   once per response, distributional parameter and
#                   nonlinear parameter, which is how brms keys its sd
#                   rows and how frmtmb resolves class "sd"
#   cor             lkj(1), uniform over correlation matrices, carried
#                   onto frmtmb's own correlation parameters with the
#                   Jacobian of that map (see R/priors.R)
#   rescor          lkj(1) on the residual correlation of
#                   set_rescor(TRUE), and on a logistic_normal()
#                   response's lncor; its sigma<k> take the sigma rows
#                   below with s = 2.5
#   sigma           student_t(3, 0, s) on the natural scale, when sigma
#                   is intercept-only; student_t(3, 0, 2.5) on the LOG
#                   scale when sigma carries its own predictor (brms
#                   makes the same distinction)
#   where s = max(2.5, round(mad(y*), 1)), and a multivariate model reads
#   y* and s off each response for that response's rows
#
# The transform is `brms:::def_scale_prior.brmsterms()`, followed line
# for line rather than approximated:
#
#   - the transformable links are identity, log, inverse, sqrt and
#     1/mu^2, and a LOG-SCALE family (lognormal and its relatives, whose
#     mu link is spelled identity but whose response is logged inside
#     the density) counts as log whatever its link says. Every other
#     link - logit, cloglog, ... - keeps location 0 and scale 2.5, which
#     is what brms reports for a bernoulli model;
#   - a response that was a FACTOR, or whose family spreads its dpars
#     over categories (ordinal, categorical, multinomial), is not
#     transformed either, so those keep 0 and 2.5. brms spells that
#     `!is_like_factor(y) && !conv_cats_dpars(x)`;
#   - the zero shift is PER ELEMENT and applies to zeros only -
#     `ifelse(y == 0, y + 0.1, y)` - and only under log/inverse/1/mu^2,
#     not under sqrt. Shifting the whole vector instead moves the
#     median: a count response with median 1 and some zeros gets
#     location log(1) = 0 from brms and log(1.1) = 0.1 from the naive
#     spelling;
#   - `round(mad(y*), 1)` raises the scale only when it is finite, and
#     `round(median(y*), 1)` sets the location only when it is finite.
#
# WHAT IS NOT MATCHED, and why, is documented in ?frm_sample and said
# out loud by the disclosure message: a correlation whose STRUCTURE has
# no LKJ density (toep(), which is not positive definite everywhere, so
# there is no correlation matrix over its parameter space to put one
# on), sum-to-zero ordinal thresholds (brms's prior sits on a vector
# frmtmb does not have), and NONLINEAR parameters (class "b" in brms,
# which brms leaves flat).

# Families whose response is on the log scale inside the density even
# though the mu link is spelled identity; brms's `has_logscale()`.
logscale_families <- c("lognormal", "shifted_lognormal",
                       "hurdle_lognormal")

#' brms's `def_scale_prior()` location and scale for one response.
#'
#' The location is the INTERCEPT's, before any offset: brms then takes
#' the mean offset of the predictor the intercept belongs to away from
#' it (`default_intercept_location()`).
#'
#' `link` is the link brms transforms the response by. brms reads the
#' family's own link for the `sd` rows and for the location, and swaps
#' in the dpar's link for a dpar: `"identity"` for one without a
#' formula, its link for one with a formula (`def_dpar_prior()`), and a
#' mixture component's own link for its location. A mixture family's
#' link is the identity, and brms does not see a log-scale component
#' through it.
#'
#' @noRd
default_prior_scale <- function(fit, rspec = fit$spec$responses[[1L]],
                                link = NULL) {
  fam <- rspec$family
  y <- fit$frame[["y"]][[rspec$resp_name]]
  mix <- !is.null(fam[["mix"]])
  link <- if (!mix && fam[["family"]] %in% logscale_families) "log" else {
    link %||% fam[["links"]][["mu"]]$name %||% "identity"
  }
  flat <- list(location = 0, scale = 2.5, link = link, centered = FALSE)
  if (is.null(y) || !is.numeric(y) || is.matrix(y) || !length(y)) {
    return(flat)
  }
  # a factor response (its codes are labels, not numbers) and a family
  # whose dpars run over categories both skip the transform
  if (!is.null(fit$frame[["y_levels"]][[rspec$resp_name]])) return(flat)
  if (fam[["type"]] %in% c("ordinal", "categorical")) return(flat)
  if (!link %in% c("identity", "log", "inverse", "sqrt", "1/mu^2")) {
    return(flat)
  }
  if (link %in% c("log", "inverse", "1/mu^2")) {
    y <- ifelse(y == 0, y + 0.1, y)
  }
  # a negative response under a log or sqrt link gives NaN, which the
  # finiteness checks below read as "keep 2.5"; brms silences the same
  # warning (SW() in def_scale_prior())
  yl <- suppressWarnings(switch(link, identity = y, log = log(y),
                                inverse = 1 / y, sqrt = sqrt(y),
                                `1/mu^2` = 1 / y^2))
  loc <- round(stats::median(yl), 1L)
  scl <- round(stats::mad(yl), 1L)
  list(location = if (is.finite(loc)) loc else 0,
       scale = if (is.finite(scl)) max(2.5, scl) else 2.5,
       link = link, centered = TRUE)
}

#' brms's default location for one predictor's intercept: the
#' response's location less the predictor's mean offset. brms's
#' `def_scale_prior()` subtracts it because an offset moves the
#' intercept by that much (`y ~ 1 + offset(o)` with `y` near 2 and `o` of
#' 10 puts the intercept near -8), and only where the response was
#' transformed at all: a logit-link model keeps location 0.
#'
#' @noRd
default_intercept_location <- function(ps, lp) {
  off <- lp[["offset"]]
  if (!isTRUE(ps$centered) || !length(off)) return(ps$location)
  m <- mean(off)
  if (is.finite(m)) ps$location - m else ps$location
}

#' A `set_prior()` spec whose distribution sits on `exp(coefficient)`
#' rather than on the coefficient, with the class `"sd"` log-Jacobian.
#'
#' @noRd
natural_dpar_prior <- function(dist, dpar, resp = "", lb = NA, ub = NA) {
  spec <- unclass(set_prior(dist, class = "Intercept", dpar = dpar,
                            resp = resp, lb = lb, ub = ub))
  spec[[1L]]$natural <- TRUE
  structure(spec, class = "frmtmb_priorlist")
}

#' The prefixes brms keys its class `"sd"` rows by, one per standard
#' deviation a block carries: response (multivariate only),
#' distributional parameter (`""` for the location ones) and nonlinear
#' parameter. frmtmb resolves a class `"sd"` specification by the same
#' prefix, so a default has to be written once per prefix to reach
#' every standard deviation, as brms writes one per row.
#'
#' brms's default length-scale prior for each `gp()` sub-GP, its
#' `def_lscale_prior()`: an inverse gamma putting 1 percent of its mass
#' below the smallest distance between two of the sub-GP's scaled
#' positions (floored at 1 percent of the largest) and 1 percent above
#' the largest, one per covariate when the term is not isotropic.
#' `normal(0, 0.5)` where the two conditions have no solution, as in
#' brms. A list of `set_prior()` arguments.
#'
#' @noRd
default_lscale_priors <- function(fit) {
  multi <- length(fit$spec$responses) > 1L
  out <- list()
  for (lp in fit$frame[["linpreds"]]) {
    rspec <- fit$spec$responses[[lp[["resp"]]]]
    dp <- lp[["dpar"]]
    for (gi in lp[["gps"]] %||% list()) {
      bk <- fit$frame[["re_blocks"]][[gi[["block_id"]]]]
      X <- gi[["prior_x"]]
      cf <- bk[["gp_brms"]][["sfx2"]]
      if (is.null(X) || is.null(cf)) next
      pr <- if (isTRUE(bk[["gp_iso"]])) {
        brms_lscale_prior(X)
      } else {
        vapply(seq_len(ncol(X)), function(j) {
          brms_lscale_prior(X[, j, drop = FALSE])
        }, "")
      }
      for (j in seq_along(cf)) {
        out[[length(out) + 1L]] <- list(
          prior = pr[j], coef = cf[j],
          resp = if (multi) lp[["resp"]] else "",
          dpar = default_lp_dpar(rspec, dp),
          nlpar = if (dp %in% (rspec$nlpars %||% character(0))) dp else "")
      }
    }
  }
  out
}

#' One `def_lscale_prior()` density, brms 2.23.0's own construction:
#' `nleqslv()` from `c(0, 0)` on the log shape and log scale. brms reads
#' the answer only when the solver converged and rounds it to six
#' decimals, which makes any solver that converges to the unique root
#' give the same string; this one is Newton's method with a central
#' difference Jacobian and step halving.
#'
#' @noRd
brms_lscale_prior <- function(X, plb = 0.01, pub = 0.01) {
  X <- as.matrix(X)
  dq <- 0
  for (j in seq_len(ncol(X))) dq <- dq + outer(X[, j], X[, j], "-")^2
  ub <- sqrt(max(dq))
  lb <- sqrt(min(dq[dq > 0]))
  lb <- max(lb, 0.01 * ub)
  # brms's pinvgamma(q, shape, scale) is pgamma(1 / q, shape, rate =
  # scale) in the opposite tail
  fn <- function(x) {
    a <- exp(x[1])
    b <- exp(x[2])
    c(stats::pgamma(1 / lb, a, rate = b, lower.tail = FALSE, log.p = TRUE) -
        log(plb),
      stats::pgamma(1 / ub, a, rate = b, lower.tail = TRUE, log.p = TRUE) -
        log(pub))
  }
  x <- c(0, 0)
  f <- fn(x)
  ok <- FALSE
  for (it in seq_len(200L)) {
    if (!all(is.finite(f))) break
    if (max(abs(f)) < 1e-12) {
      ok <- TRUE
      break
    }
    J <- vapply(1:2, function(k) {
      h <- c(0, 0)
      h[k] <- 1e-6
      (fn(x + h) - fn(x - h)) / 2e-6
    }, numeric(2))
    step <- tryCatch(solve(J, -f), error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) break
    t <- 1
    repeat {
      xn <- x + t * step
      fnew <- fn(xn)
      if (all(is.finite(fnew)) && sum(fnew^2) < sum(f^2)) break
      t <- t / 2
      if (t < 1e-10) break
    }
    if (t < 1e-10) break
    x <- xn
    f <- fnew
  }
  if (!ok) return("normal(0, 0.5)")
  p <- round(exp(x), 6)
  sprintf("inv_gamma(%s, %s)", format(p[1], digits = 15),
          format(p[2], digits = 15))
}

#' `gp = TRUE` gives the same prefixes for the `gp()` blocks, which brms
#' priors as class `"sdgp"` rather than `"sd"`.
#'
#' @noRd
default_sd_prefixes <- function(fit, gp = FALSE) {
  multi <- length(fit$spec$responses) > 1L
  out <- list()
  for (bk in fit$frame[["re_blocks"]]) {
    if (!length(block_sd_idx(bk))) next
    if (xor(gp, bk[["covstruct"]] %in% c("gp", "hsgp"))) next
    for (cp in bk[["components"]] %||% list()) {
      lp <- fit$frame[["linpreds"]][[cp[["lp_key"]] %||% ""]]
      if (is.null(lp)) next
      rspec <- fit$spec$responses[[lp[["resp"]]]]
      nl <- rspec$nlpars %||% character(0)
      dp <- lp[["dpar"]]
      key <- c(resp = if (multi) lp[["resp"]] else "",
               dpar = default_lp_dpar(rspec, dp),
               nlpar = if (dp %in% nl) dp else "")
      out[[paste(key, collapse = "|")]] <- c(key, rname = lp[["resp"]])
    }
  }
  unname(out)
}

#' The `dpar` a default row carries: empty for the location of a family
#' that has one and for a nonlinear parameter, the parameter's name
#' otherwise. A categorical, multinomial or mixture family has several
#' locations, and brms names each (`mub`, `mu1`). frmtmb names them for
#' every other family with several location dpars too, and refuses a
#' location prior without that name (`lp_prior_dpar()`, which this
#' mirrors). Before, one specification per location was written without
#' it, and each reached every location.
#'
#' @noRd
default_lp_dpar <- function(rspec, dp) {
  if (dp %in% (rspec$nlpars %||% character(0))) return("")
  loc <- setdiff(rspec$primary_dpars, rspec$nlpars %||% character(0))
  several <- length(loc) > 1L || !is.null(rspec$family[["mix"]]) ||
    isTRUE(rspec$family[["family"]] %in% c("categorical", draws_count_cats))
  if (dp %in% rspec$primary_dpars && !several) "" else dp
}

#' Whether `lp` is a `gen_extreme_value()` `xi` without a formula. With
#' a formula, `xi ~ 1` included, brms samples `xi` on its `log1p` link
#' with `normal(0, 4)` on the intercept, which `brms_def_dpar_prior()`
#' gives it. The two spellings have the same design, one intercept
#' column, so the tell is the formula, which `dpar_has_predictor()`
#' reads; a design test gave `xi ~ 1` the scalar's note and no default
#' (review fix 6).
#'
#' @noRd
gev_xi_scalar <- function(spec, rspec, lp) {
  identical(rspec$family[["family"]], "gen_extreme_value") &&
    identical(lp[["dpar"]], "xi") && is.null(lp[["constant"]]) &&
    is.null(lp[["equate"]]) && !dpar_has_predictor(spec, lp)
}

#' brms's class for one of a response's dpars, its `dpar_class()`: a
#' mixture component's `sigma2` is class `sigma` and its `theta2` class
#' `theta`, and a logistic_normal() category's `sigma<k>` is class
#' `sigma`. frmtmb calls com_poisson()'s dispersion `nu` where brms calls
#' it `shape`, and brms's default belongs to the quantity, not the word.
#'
#' @noRd
brms_dpar_class <- function(rspec, dp) {
  fam <- rspec$family
  if (dp %in% fam[["simplex"]][["sigma"]]) return("sigma")
  if (!is.null(fam[["mix"]])) return(sub("[0-9]+$", "", dp))
  if (identical(fam[["family"]], "compois") && identical(dp, "nu")) {
    return("shape")
  }
  dp
}

#' The defaults a brms family states for itself (brms 2.23.0's
#' `family$prior` hooks, brms#1614), or `NULL` where the family has none
#' for the class under that link. A mixture never reaches them: brms
#' asks the mixture family, which has none.
#'
#' @noRd
brms_family_prior <- function(rspec, cls, link) {
  fam <- rspec$family
  if (!is.null(fam[["mix"]])) return(NULL)
  switch(fam[["family"]],
    xbeta = if (identical(cls, "kappa")) {
      if (identical(link, "identity")) "gamma(0.01, 0.01)" else
        "student_t(3, 0, 2.5)"
    },
    negbinomial = , hurdle_negbinomial = , zero_inflated_negbinomial =
      if (identical(cls, "shape") && identical(link, "identity")) {
        "inv_gamma(0.4, 0.3)"
      },
    # the same prior moved to sigma = 1 / shape: brms's hook for
    # negbinomial2 says gamma where negbinomial's says inv_gamma
    negbinomial2 = if (identical(cls, "sigma") && identical(link, "identity")) {
      "gamma(0.4, 0.3)"
    },
    von_mises = if (identical(cls, "mu") && identical(link, "tan_half")) {
      "student_t(1, 0, 1)"
    },
    NULL)
}

#' brms 2.23.0's `def_dpar_prior()` for one dpar of class `cls`, as a
#' density string, `""` where brms leaves the slot flat. `link` is
#' `"identity"` for a dpar without a formula, whose density is then on
#' the parameter itself, and the dpar's link for one with a formula,
#' whose density is then on its intercept. A family brms does not have
#' gets the class's default, which is what brms gives a custom family
#' with a dpar of that name. `lk` is the dpar's link object, which the
#' `ndt` default reads (`ndt_upper()`).
#'
#' @noRd
brms_def_dpar_prior <- function(fit, rspec, cls, link, lk = NULL) {
  hook <- brms_family_prior(rspec, cls, link)
  if (!is.null(hook)) return(hook)
  if (identical(cls, "sigma")) {
    ps <- default_prior_scale(fit, rspec, link = link)
    return(sprintf("student_t(3, 0, %s)", as.character(ps$scale)))
  }
  if (identical(link, "identity")) {
    switch(cls,
      shape = "gamma(0.01, 0.01)", nu = "gamma(2, 0.1)",
      phi = "gamma(0.01, 0.01)", kappa = "gamma(2, 0.01)",
      beta = "gamma(1, 0.1)", zi = , hu = , zoi = , coi = , bias = ,
      quantile = "beta(1, 1)", bs = "gamma(1, 1)",
      ndt = ndt_uniform_prior(fit, rspec, lk), xi = "normal(0, 2.5)",
      alpha = "normal(0, 4)", disc = "lognormal(0, 1)",
      theta = "logistic(0, 1)", "")
  } else {
    switch(cls,
      shape = "student_t(3, 0, 2.5)", nu = "normal(2.7, 0.8)",
      phi = "student_t(3, 0, 2.5)", kappa = "normal(5.0, 0.8)",
      beta = "normal(1.7, 1.3)", zi = , hu = , zoi = , coi = , bias = ,
      quantile = "logistic(0, 1)", bs = "normal(-0.6, 1.3)",
      xi = "normal(0, 4)", alpha = "normal(0, 4)", disc = "normal(0, 1)",
      "")
  }
}

#' brms's bounds on a distributional parameter's own class
#' (`dpar_bounds()`), where frmtmb's link does not hold the parameter
#' inside them already: `nu` above 1 under an identity or log link, and
#' a log-linked `ndt` below the smallest response. brms declares the
#' bound on the parameter (`real<lower=1> nu`,
#' `real<lower=0, upper=min_Y> ndt`), so Stan's transform keeps a draw
#' inside it; without it the density is a wall the sampler meets
#' (dev/reviews/2026-10-14-sampleprior.md, findings 3 and 4). A link
#' whose image is inside the bound, logit or logm1 or frmtmb.eam's
#' bounded `ndt`, needs none and gets none.
#'
#' @noRd
brms_natural_bounds <- function(fit, rspec, cls, link) {
  bd <- switch(cls,
    sigma = , shape = , phi = , kappa = , beta = , bs = , disc = c(0, NA),
    nu = c(1, NA),
    zi = , hu = , zoi = , coi = , bias = , quantile = c(0, 1),
    ndt = c(0, ndt_upper(fit, rspec, link)),
    c(NA, NA))
  img <- tryCatch(suppressWarnings(link$linkinv(c(-30, 30))),
                  error = function(e) c(-Inf, Inf))
  lb <- if (!is.na(bd[1L]) && isTRUE(img[1L] < bd[1L])) bd[1L] else NA
  ub <- if (!is.na(bd[2L]) && isTRUE(img[2L] > bd[2L])) bd[2L] else NA
  c(lb = lb, ub = ub)
}

#' The smallest response, brms's `min_Y`, or `NA` without a positive one.
#'
#' @noRd
ndt_min_y <- function(fit, rspec) {
  y <- fit$frame[["y"]][[rspec$resp_name]]
  if (is.matrix(y)) y <- y[, 1L]
  y <- y[is.finite(y)]
  if (!length(y) || !(min(y) > 0)) return(NA_real_)
  min(y)
}

#' The top of a non-decision time's support: the top of its link's
#' image where the link is bounded above, and brms's `min_Y`, the
#' smallest response, where it is not. frmtmb.eam bounds `ndt` in its
#' link. Without `ndt_group()` that bound is the fastest response, so
#' the image is `(0, min_Y)`, brms's interval; with `ndt_group()` the
#' parameter is a FRACTION of each group's fastest response, the image
#' is `(0, 1)`, and `min_Y` would cap the fraction where it does not
#' belong (dev/reviews/2026-10-14-sampleprior.md, finding 9: the
#' posterior pinned at 0.33 against an ML fraction of 0.818). A log,
#' softplus or identity link, the shifted lognormal's, is unbounded and
#' takes `min_Y`.
#'
#' @noRd
ndt_upper <- function(fit, rspec, link) {
  top <- if (is.null(link)) Inf else {
    tryCatch(suppressWarnings(link$linkinv(Inf)), error = function(e) NaN)
  }
  if (!is.finite(top)) {
    # a linkinv() that cannot take Inf: read it far out instead
    top <- tryCatch(suppressWarnings(link$linkinv(800)),
                    error = function(e) Inf)
  }
  if (is.finite(top) && top > 0) top else ndt_min_y(fit, rspec)
}

#' brms's `uniform(0, min_Y)` on a non-decision time without a formula,
#' with `min_Y` the smallest response, over the link's own image where
#' the link is bounded above (`ndt_upper()`). The bound is written to 17
#' significant digits, so the density's edge is the support's own value
#' and not a rounding of it.
#'
#' @noRd
ndt_uniform_prior <- function(fit, rspec, link = NULL) {
  m <- ndt_upper(fit, rspec, link)
  if (is.na(m)) return("")
  sprintf("uniform(0, %.17g)", m)
}

#' brms's default priors for the model this object holds, as a
#' `frmtmb_priorlist`, or `NULL` when the model has no slot they cover.
#'
#' A multivariate model takes each response's defaults from that
#' response, with `resp` set, as brms's rows do; `cor` and `rescor` are
#' class-wide.
#'
#' @noRd
default_priors_for <- function(fit) {
  multi <- length(fit$spec$responses) > 1L
  pl <- NULL
  add <- function(p) pl <<- if (is.null(pl)) p else pl + p
  st <- function(loc, scl) {
    sprintf("student_t(3, %s, %s)", as.character(loc), as.character(scl))
  }
  scales <- lapply(fit$spec$responses, function(rs) {
    default_prior_scale(fit, rs)
  })

  for (lp in fit$frame[["linpreds"]]) {
    # an equated dpar (sigma1 = "sigma2") owns no parameter: brms puts
    # the prior on its target only
    if (!is.null(lp[["constant"]]) || !is.null(lp[["nl_body"]]) ||
          !is.null(lp[["equate"]])) next
    if (!"(Intercept)" %in% colnames(lp[["X"]])) next
    # `0 + Intercept` or center = FALSE: brms's intercept is then class
    # "b", which brms leaves flat
    if (isFALSE(lp[["center"]])) next
    rspec <- fit$spec$responses[[lp[["resp"]]]]
    # the shared nu of a Student-t rescor model has no response, as its
    # brms row has none
    rs <- if (multi && !isTRUE(lp[["shared"]])) lp[["resp"]] else ""
    dp <- lp[["dpar"]]
    # a NONLINEAR parameter's coefficients are class "b" in brms, and
    # brms leaves class "b" flat: the response's median and mad say
    # nothing about a rate or a shape sitting inside a nonlinear body,
    # so a default there would be an invented prior rather than brms's
    if (dp %in% (rspec$nlpars %||% character(0))) next
    if (dp %in% rspec$primary_dpars) {
      # a mixture component's location is read on its own link
      ps <- if (is.null(rspec$family[["mix"]])) scales[[lp[["resp"]]]] else {
        default_prior_scale(fit, rspec, link = lp[["link"]]$name)
      }
      d <- brms_family_prior(rspec, "mu", lp[["link"]]$name) %||%
        st(default_intercept_location(ps, lp), ps$scale)
      add(set_prior(d, class = "Intercept", resp = rs,
                    dpar = default_lp_dpar(rspec, dp)))
      next
    }
    if (gev_xi_scalar(fit$spec, rspec, lp)) {
      # brms's default on a scalar xi is normal(0, 2.5) on tmp_xi, and
      # frm_sample() samples that xi through tmp_xi as brms does
      # (gev_xi_chart_frame()), where this density is read on tmp_xi.
      # Where the chart cannot be used it is a density on xi itself
      add(set_prior("normal(0, 2.5)", class = "xi", resp = rs))
      next
    }
    # every other dpar: brms's def_dpar_prior(), which reads the dpar's
    # link only when it has a formula (`dp ~ 1` included); without one
    # the density is on the parameter itself, in the dpar's own class
    pred <- dpar_has_predictor(fit$spec, lp)
    cls <- brms_dpar_class(rspec, dp)
    # a mixture's weights without a formula take brms's dirichlet(1) on
    # the whole simplex instead (mix_theta_nlp())
    if (identical(cls, "theta") && !pred) next
    d <- brms_def_dpar_prior(fit, rspec, cls,
                             if (pred) lp[["link"]]$name else "identity",
                             lk = lp[["link"]])
    if (!nzchar(d)) next
    add(if (pred) {
      set_prior(d, class = "Intercept", dpar = dp, resp = rs)
    } else {
      bd <- brms_natural_bounds(fit, rspec, cls, lp[["link"]])
      p <- natural_dpar_prior(d, dp, resp = rs, lb = bd[["lb"]],
                              ub = bd[["ub"]])
      # brms's sigma under a first-order cov = TRUE ar(), ma() or arma()
      # is the innovation sd, and core reads every density on class sigma
      # there (autocor_innov()); the flag only marks the message line
      if (identical(cls, "sigma") &&
            lp[["resp"]] %in% autocor_innov_sd(fit)) {
        p[[1L]]$innov_sd <- TRUE
      }
      p
    })
  }
  # an ordinal model's thresholds are brms's class "Intercept", with the
  # response's def_scale_prior(): one row per threshold vector, keyed as
  # get_prior() keys it, and none under sum_to_zero, where brms's
  # density is on a vector frmtmb does not have
  for (ob in ord_lp_blocks(fit$frame, fit$spec)) {
    if (ob$shared && !identical(ob$k, 1L)) next
    if (!length(fit$frame[["par_template"]][[ob$comp]] %||% numeric(0))) {
      next
    }
    if (identical(ob$fam[["thres"]][["type"]], "sum_to_zero")) next
    r <- ob$lp[["resp"]] %||% names(fit$spec$responses)[1L]
    add(set_prior(st(0, scales[[r]]$scale), class = "Intercept",
                  dpar = ob$prior_dpar, resp = if (multi) r else ""))
  }

  for (k in default_sd_prefixes(fit)) {
    add(set_prior(st(0, scales[[k[["rname"]]]]$scale), class = "sd",
                  resp = k[["resp"]], dpar = k[["dpar"]],
                  nlpar = k[["nlpar"]]))
  }
  # a gp() term's sd: brms's prior_gp() writes the same def_scale_prior
  # on class "sdgp", and its length scales an inv_gamma tuned to each
  # sub-GP's own distances (def_lscale_prior())
  for (k in default_sd_prefixes(fit, gp = TRUE)) {
    add(set_prior(st(0, scales[[k[["rname"]]]]$scale), class = "sdgp",
                  resp = k[["resp"]], dpar = k[["dpar"]],
                  nlpar = k[["nlpar"]]))
  }
  for (e in default_lscale_priors(fit)) {
    add(set_prior(e$prior, class = "lscale", coef = e$coef, resp = e$resp,
                  dpar = e$dpar, nlpar = e$nlpar))
  }
  # lkj(1), brms's own default: uniform over correlation matrices. Only
  # when some block's correlation HAS an LKJ density, so that a model
  # whose only correlated block is a toep() does not fail on a default
  # it cannot honor (default_prior_notes() names that block instead)
  has_cor <- any(vapply(fit$frame[["re_blocks"]], function(bk) {
    identical(block_cor_prior(bk), "lkj")
  }, TRUE))
  if (has_cor) add(set_prior("lkj(1)", class = "cor"))
  # brms's lkj(1) on the residual correlation of set_rescor(TRUE)
  if (length(fit$frame[["par_template"]][["thetar"]] %||% numeric(0))) {
    add(set_prior("lkj(1)", class = "rescor"))
  }
  # and on a logistic_normal() response's correlations, class lncor,
  # one per such response; a multivariate model holds each under its
  # own name (core's extra_map) and names its response, as brms does
  for (r in names(fit$spec$responses)) {
    if (!identical(fit$spec$responses[[r]]$family[["simplex"]]$kind,
                   "logistic_normal")) next
    comp <- fit$frame[["extra_map"]][[r]][["lncor_raw"]] %||% "lncor_raw"
    if (!length(fit$frame[["par_template"]][[comp]] %||% numeric(0))) next
    add(set_prior("lkj(1)", class = "lncor", resp = if (multi) r else ""))
  }
  # brms's (flat) on a residual-correlation parameter frmtmb holds
  # through a transform (a cov = TRUE ar or ma coefficient, cosy, sar()'s
  # lagsar): a uniform density on the natural parameter over the whole
  # image of the map, which core carries inward with the map's
  # log-Jacobian, so rho is flat where it lives and not on the scale the
  # term parameterizes it on, which for ar, ma and cosy was the default
  # once: (1 - rho^2)^(-3/2) on a first-order rho, with its mass at the
  # ends. cosy's lb = 0 is the bound brms declares inside its window
  ap <- autocor_flat_priors(fit)
  for (i in seq_len(nrow(ap))) {
    p <- set_prior(sprintf("uniform(%.17g, %.17g)", ap$lower[i],
                           ap$upper[i]), class = ap$class[i],
                   resp = if (multi) ap$resp[i] else "", lb = ap$lb[i],
                   ub = ap$ub[i])
    p[[1L]]$flat_interval <- TRUE
    # a higher-order cov = TRUE ar or ma, which brms does not have: the
    # uniform is on a box around the stationary region, not an interval
    p[[1L]]$flat_region <- !ap$interval[i]
    add(p)
  }
  # and on an unstr() residual correlation, class cortime, per response
  for (r in names(fit$frame[["autocor"]] %||% list())) {
    ac <- fit$frame[["autocor"]][[r]]
    if (identical(ac[["struct"]], "unstr") && isTRUE(ac[["d"]] >= 2L)) {
      add(set_prior("lkj(1)", class = "cortime", resp = if (multi) r else ""))
    }
  }
  pl
}

#' The mixing weights of each mixture response that has no formula on
#' any of them, where brms 2.23.0 puts `dirichlet(1)` on the simplex
#' (`has_joint_theta()`): for each, the positions of its `K - 1` log
#' ratios in `betad`. A weight with a formula is brms's class
#' `Intercept` instead, and gets `logistic(0, 1)` there.
#'
#' @noRd
mix_theta_terms <- function(fit) {
  tab <- brms_coef_table(fit)
  smp <- attr(tab, "simplex")
  if (!length(smp)) return(list())
  n_beta <- length(fit$frame[["par_template"]][["beta"]])
  n_bd <- length(fit$frame[["par_template"]][["betad"]])
  free <- setdiff(seq_len(n_bd), fit$frame[["betad_fixed_idx"]])
  lapply(smp, function(s) {
    list(idx = free[s$pos - n_beta], K = length(s$names),
         names = s$names)
  })
}

#' brms's `dirichlet(1)` on a mixture's weights, as a negative log
#' density of the log ratios the fit holds them in. The flat density on
#' the simplex carried to those coordinates is the softmax's Jacobian,
#' `Gamma(K) prod(w)` over all `K` weights, whichever component is the
#' reference. `NULL` without such weights.
#'
#' @noRd
mix_theta_nlp <- function(fit) {
  terms <- if (isTRUE(attr(fit$prior, "mix_theta"))) mix_theta_terms(fit)
  if (!length(terms)) return(NULL)
  function(pars) {
    "c" <- RTMB::ADoverload("c")
    out <- 0
    for (tm in terms) {
      a <- c(0, pars[["betad"]][tm$idx])
      lse <- a[1L]
      for (k in seq_along(a)[-1L]) lse <- RTMB::logspace_add(lse, a[k])
      out <- out - (sum(a) - tm$K * lse) - lgamma(tm$K)
    }
    out
  }
}

#' What this model has that brms would prior and frmtmb deliberately
#' leaves flat. Every such slot is named out loud in the disclosure
#' message, so a family that gets no defaults never passes in silence.
#' `user_pl` is the prior the run applies besides the defaults: a
#' simplex with a Dirichlet row of its own there does not get brms's
#' default, so the note leaves it out.
#'
#' @noRd
default_prior_notes <- function(fit, user_pl = NULL) {
  multi <- length(fit$spec$responses) > 1L
  notes <- character(0)
  for (rspec in fit$spec$responses) {
    # a multivariate model says which response each note is about
    of <- if (multi) paste0(" of response ", rspec$resp_name) else ""
    nlp <- rspec$nlpars %||% character(0)
    if (length(nlp)) {
      notes <- c(notes, paste0("nonlinear parameters", of, " stay flat: ",
                               paste(nlp, collapse = ", "),
                               " (brms leaves them flat too, as class b; ",
                               "set_prior(nlpar = ) writes them)"))
    }
    for (ob in ord_lp_blocks(fit$frame, fit$spec)) {
      r <- ob$lp[["resp"]] %||% rspec$resp_name
      if (!identical(r, rspec$resp_name)) next
      if (!identical(ob$fam[["thres"]][["type"]], "sum_to_zero")) next
      notes <- c(notes, paste0(
        "no default for the sum_to_zero thresholds", of, " (brms puts ",
        "student_t(3, 0, 2.5) on thresholds it declares before centering ",
        "them, which frmtmb does not have)"))
      break
    }
    for (tm in gev_xi_chart_terms(fit)) {
      if (!identical(tm$resp, rspec$resp_name)) next
      notes <- c(notes, if (length(tm$blocked)) {
        gev_xi_blocked_text(tm, of)
      } else {
        paste0("xi", of, " is sampled as brms samples a scalar xi, ",
               "through tmp_xi, mapped into the support each draw's mu and ",
               "sigma leave, and above -1, where brms's range goes on; a ",
               "prior on class xi is a prior on tmp_xi, as in brms")
      })
    }
  }
  # a mixture's weights without a formula: brms's dirichlet(1) on the
  # simplex, which no set_prior() class addresses here
  for (tm in mix_theta_terms(fit)) {
    notes <- c(notes, paste0(
      "theta: dirichlet(1), brms's default, on the mixing weights (",
      paste(tm$names, collapse = ", "), "); it is the flat density on ",
      "the simplex, and prior = \"flat\" drops it with the rest"))
  }
  # correlations get lkj(1), brms's own default, wherever a density
  # exists for the block's parameterization; what is left is named
  ungated <- unique(vapply(fit$frame[["re_blocks"]], function(bk) {
    if (identical(block_cor_prior(bk), "unsupported")) {
      paste0(bk[["term_label"]], " [", bk[["covstruct"]], "]")
    } else {
      ""
    }
  }, ""))
  ungated <- ungated[nzchar(ungated)]
  if (length(ungated)) {
    notes <- c(notes, paste0("no correlation defaults for ",
                             paste(ungated, collapse = ", "),
                             " (no LKJ density fits its parameters)"))
  }
  # the simplexes the user's own rows prior take that density instead
  # (mo_simplex_nlp()), so naming them here misreported the run
  own <- if (!is.null(user_pl)) {
    ent <- tryCatch(resolve_prior_input(fit, user_pl)$entries,
                    error = function(e) list())
    unlist(lapply(ent, function(e) {
      if (identical(e$dist$kind, "dirichlet")) e$comp
    }))
  }
  mo <- Filter(function(t) !t$zeta %in% own, mo_frame_terms(fit))
  if (length(mo)) {
    notes <- c(notes, paste0(
      "simo: dirichlet(1), brms's default, on each mo() simplex without ",
      "a simo row of its own (",
      paste(sub("[[]1[]]$", "", vapply(mo, function(t) t$names[1L], "")),
            collapse = ", "),
      "); it is the flat density on the simplex, and prior = \"flat\" ",
      "keeps it"))
  }
  # every simplex a family declares (frmtmb_family(simplex_pars =)), the
  # cox() baseline's sbhaz among them, by the brms name of its draws; a
  # mo() term carries no `kind`, and the note above names it
  fs <- Filter(function(t) !is.null(t$kind) && !t$zeta %in% own,
               simplex_frame_terms(fit))
  seen <- character()
  for (t in fs) {
    stem <- sub("[[]1[]]$", "", t$names[1L])
    if (stem %in% seen) next
    seen <- c(seen, stem)
    notes <- c(notes, paste0(
      stem, ": dirichlet(1), brms's default, on the simplex the ", t$kind,
      "() family declares; it is the flat density on the simplex, and ",
      "prior = \"flat\" keeps it"))
  }
  notes
}

#' The user specifications whose every parameter keeps a default, one
#' note each for the announcement.
#'
#' A default steps aside only for its own slot (`drop_superseded()`),
#' and a default written per coefficient is more specific than a
#' class-wide user row, so it wins on its parameter, as brms's does:
#' `set_prior("inv_gamma(2, 1)", class = "lscale")` on a `gp(x)` model
#' keeps brms's tuned `inv_gamma()` on `lscale_gpx`
#' (dev/prilabel-rev-lscale.R). A user's own coef row shadows a
#' class-wide row of theirs the same way. Such a specification applies
#' nowhere, which would pass in silence; the announcement names it and
#' the coef that reaches the parameter.
#'
#' @noRd
unreached_prior_notes <- function(fit, user_pl, ri) {
  if (is.null(user_pl) || is.null(ri)) return(character(0))
  key <- function(e) paste0(e$comp, ".", e$idx)
  got <- list()
  for (e in ri$entries) {
    for (k in key(e)) got[[k]] <- e$dist
  }
  out <- character(0)
  for (s in unclass(user_pl)) {
    if (is.null(s$dist)) next
    one <- structure(list(s), class = "frmtmb_priorlist")
    alone <- tryCatch(resolve_prior_input(fit, one)$entries,
                      error = function(e) NULL)
    if (!length(alone)) next
    # the density the specification alone puts at a position, as the
    # resolver writes it (an lkj or a Dirichlet is rewritten per block),
    # against what the whole stack leaves there
    hit <- vapply(alone, function(e) {
      any(vapply(key(e), function(k) identical(got[[k]], e$dist), NA))
    }, NA)
    if (any(hit)) next
    out <- c(out, paste0(
      spec_target(s), " applies to no parameter: a more specific prior ",
      "(a default or another specification) on each of them stays, as ",
      "in brms. To replace a default there, name the parameter's coef ",
      "as get_prior(route = \"sample\") lists it"))
  }
  out
}

#' Drop the specifications a later one for the same slot replaces:
#' same class, coef, group, response, dpar and nlpar. This is brms's
#' `update(prior = )`, which puts a new row in place of the stored row
#' for its slot, and it is what keeps `prior_summary()` from listing
#' both a MAP fit's own prior and a `prior =` on the sampling call.
#'
#' A bounds-only specification replaces the stored one too, density and
#' all, as the brms row it mirrors does; the call writes the whole slot.
#'
#' @noRd
drop_superseded <- function(base, over) {
  tg <- vapply(unclass(over), spec_target, "")
  keep <- Filter(function(s) !(spec_target(s) %in% tg), unclass(base))
  if (!length(keep)) NULL else structure(keep,
                                         class = "frmtmb_priorlist")
}

#' What a user's specification of a slot keeps from brms's default
#' there, since brms keeps both. Every default with a bound keeps that
#' bound on each side the specification leaves unset: brms declares
#' `real<lower=1> nu`, `real<lower=0> sigma` under an identity link,
#' `real<lower=0,upper=min_Y> ndt` and `real<lower=0,upper=1> cosy`,
#' and a user prior changes the density, not the declaration. Through
#' frmtmb.sample 0.22.0 a user prior dropped the bound: `normal(2, 2)`
#' on nu put 3.6 percent of the draws below 1
#' (dev/review-acprior-nu.R). A flat residual-correlation default
#' (`flat_interval`) also lends a bounds-only specification its flat
#' density: brms's `set_prior("", class = "ar", lb = 0)` is flat on
#' `[0, 1]`, where without it the slot would be flat on frmtmb's
#' internal parameter.
#'
#' @noRd
keep_default_bounds <- function(user, defs) {
  flat <- Filter(function(d) {
    isTRUE(d$flat_interval) || !is.na(d$lb %||% NA) || !is.na(d$ub %||% NA)
  }, unclass(defs))
  if (!length(flat)) return(user)
  tg <- vapply(flat, spec_target, "")
  out <- lapply(unclass(user), function(s) {
    d <- flat[match(spec_target(s), tg)]
    if (!length(d) || is.null(d[[1L]])) return(s)
    d <- d[[1L]]
    if (is.na(s$lb %||% NA)) s$lb <- d$lb %||% NA
    if (is.na(s$ub %||% NA)) s$ub <- d$ub %||% NA
    if (is.null(s$dist) && isTRUE(d$flat_interval)) {
      s$dist <- d$dist
      s$prior <- d$prior
      s$flat_interval <- TRUE
      s$flat_region <- d$flat_region
    }
    s
  })
  structure(out, class = class(user))
}


#' One compact line per prior class, so the call discloses what it
#' chose - and one line per slot it deliberately left flat, so a model
#' that gets few defaults, or none, says so instead of passing in
#' silence. `message()` rather than `warning()`: this is information,
#' not a fault, and `suppressMessages()` turns it off.
#'
#' @noRd
announce_default_priors <- function(pl, notes) {
  msg <- paste0("frm_sample(): default priors (brms 2.23 defaults; ",
                "prior = \"flat\" opts out)")
  for (s in unclass(pl %||% list())) {
    kind <- if (identical(s$dist$kind, "t")) "student_t" else s$dist$kind
    a <- unlist(s$dist[-1L])
    # the ends of an interval read off eigenvalues carry rounding
    # (-0.999999999999998); the line shows them to 6 digits, and the
    # prior keeps them exact
    if (isTRUE(s$flat_interval)) {
      a <- vapply(a, function(v) format(v, digits = 6), "")
    }
    d <- paste0(kind, "(", paste(a, collapse = ", "), ")")
    # a bound brms declares inside the density's interval (cosy's 0)
    # narrows where the flat density is, so it is part of the default
    if (isTRUE(s$flat_interval)) {
      bd <- c(lb = s$lb %||% NA, ub = s$ub %||% NA)
      bd <- bd[!is.na(bd)]
      if (length(bd)) {
        d <- paste0(d, ", ", paste0(names(bd), " = ", bd, collapse = ", "))
      }
    }
    # the spelling the default is WRITTEN with, which for a density on
    # a distributional parameter itself is that parameter's own class.
    # get_prior(route = "sample") lists that word and set_prior() takes
    # it; announcing the storage pair instead named a spelling the same
    # model refuses
    sp <- spec_spelling(s)
    # and with every field that narrows it, so the line names a slot
    # set_prior() reaches: a multivariate model refuses a row without
    # its resp, and a categorical or mixture one a location without its
    # dpar
    f <- c(coef = s$coef %||% "", dpar = sp$dpar, nlpar = s$nlpar %||% "",
           resp = s$resp %||% "")
    f <- f[nzchar(f)]
    lab <- if (length(f)) {
      paste0(sp$class, " (", paste0(names(f), " = ", f, collapse = ", "),
             ")")
    } else {
      sp$class
    }
    msg <- c(msg, sprintf("  %-18s %s%s", lab, d,
                          # a scalar xi's scale has a note of its own
                          if (identical(sp$class, "xi")) ""
                          else if (isTRUE(s$flat_region))
                            "  [flat on the stationary or invertible region]"
                          else if (isTRUE(s$flat_interval))
                            "  [flat on its interval, brms's (flat)]"
                          else if (isTRUE(s$innov_sd))
                            "  [on brms's innovation sd]"
                          else if (isTRUE(s$natural)) "  [natural scale]"
                          else if (s$class %in% c("sd", "sdgp"))
                            "  [natural sd scale]"
                          else if (identical(s$class, "lscale"))
                            "  [brms's scaled inputs]"
                          else if (identical(s$class, "cor"))
                            "  [correlation matrix]" else ""))
  }
  msg <- c(msg, "  b                  (flat), as brms leaves slopes")
  for (nt in notes) {
    msg <- c(msg, paste0("  ", nt, " - see ?frm_sample"))
  }
  frm_message(paste(msg, collapse = "\n"))
  invisible(NULL)
}

#' Merge resolved prior inputs; a later one wins per parameter. Used to
#' let a user's `prior =` take over individual parameters from the
#' defaults without discarding the rest.
#'
#' @noRd
merge_prior_inputs <- function(base, over) {
  if (is.null(base)) return(over)
  if (is.null(over)) return(base)
  # collapsed, because a class "cor" entry covers a whole correlation
  # segment and its key names every position it holds
  key <- function(e) paste0(e$comp, ".", paste(e$idx, collapse = ","))
  ent <- base$entries
  names(ent) <- vapply(ent, key, "")
  for (e in over$entries) {
    # a base entry that shares ANY position with the incoming one is
    # retired, not just the one with the same key: the default LKJ term
    # spans a whole correlation segment, and leaving it in place beside
    # a per-parameter override would apply both densities there
    hit <- vapply(ent, function(b) {
      identical(b$comp, e$comp) && length(intersect(b$idx, e$idx)) > 0L
    }, TRUE)
    if (any(hit)) ent <- ent[!hit]
    ent[[key(e)]] <- e
  }
  list(entries = unname(ent),
       lower = utils::modifyList(as.list(base$lower),
                                 as.list(over$lower)),
       upper = utils::modifyList(as.list(base$upper),
                                 as.list(over$upper)))
}

#' The slots a prior specification names, for the message that has to
#' say which of the fit's own priors an opt-out discards. Both
#' spellings: a priorlist names classes, the legacy list names internal
#' parameters.
#'
#' @noRd
prior_target_list <- function(pl) {
  tg <- if (inherits(pl, "frmtmb_priorlist")) {
    vapply(unclass(pl), spec_target, "")
  } else {
    names(pl) %||% character(0)
  }
  paste(unique(tg), collapse = "; ")
}

#' Stack prior specifications least explicit first. `resolve_priorlist()`
#' already resolves a later specification over an earlier one, position
#' by position, so the stacking IS the precedence rule and nothing has
#' to be dropped by hand.
#'
#' The legacy named-list spelling addresses raw internal parameters, not
#' classes, so it cannot be concatenated into a priorlist; it rides
#' alongside in its own slot and is merged at the resolved-entry level.
#'
#' @noRd
prior_stack <- function(...) {
  pl <- NULL
  legacy <- list()
  for (p in list(...)) {
    p <- as_priorlist(p)
    if (is.null(p) || !length(p)) next
    if (inherits(p, "frmtmb_priorlist")) {
      if (!is.null(pl)) pl <- drop_superseded(pl, p)
      pl <- if (is.null(pl)) p else pl + p
    } else {
      legacy[[length(legacy) + 1L]] <- p
    }
  }
  list(pl = pl, legacy = legacy)
}

#' Resolve the priors of a sampling call: the brms defaults, the fit's
#' own prior, a specification given on the call, and the opt-out.
#'
#' PRECEDENCE, most explicit wins per addressed slot: `prior =` on the
#' call, then `base` (a fitted object's own prior, which a MAP fit
#' carries into sampling), then the defaults. The first two are stacked
#' in that order and resolved later-wins; a default steps aside only for
#' a specification of its own slot (`drop_superseded()`), and a narrower
#' one, a group's or a coefficient's, applies over it on its own
#' parameters and leaves the default on the rest. That is brms's rule:
#' `set_prior(class = "sd", group = "g", coef = "x")` keeps brms's
#' `student_t()` default on `sd_g__Intercept` and on every other group
#' (dev/prilabel-log/brms-ref.txt), where this once dropped the default
#' of the whole class.
#'
#' `defaults = FALSE` is the DIAGNOSTIC mode [check_laplace()] runs in:
#' no defaults, no disclosure, and no opt-out warning. That caller
#' withholds `base` as well, and for the same reason: a MAP fit's
#' penalty is already taped into `fit$obj`, so resolving it here would
#' rebuild a tape the caller already has, and adding anything else
#' would change the density it is measuring.
#'
#' @noRd
sample_resolve_priors <- function(fit, prior, base = NULL,
                                  defaults = TRUE) {
  if (is.character(prior)) {
    if (!identical(prior, "flat")) {
      frm_stop("prior = must be a set_prior() specification, a named list ",
               "of prior objects, or the string \"flat\" to sample the ",
               "likelihood with improper flat priors; got \"",
               paste(prior, collapse = "\", \""), "\"", call. = FALSE)
    }
    if (defaults && length(fit$frame[["re_blocks"]])) {
      frm_warning("prior = \"flat\": every variance component has a flat ",
                  "prior on its log standard deviation, under which the ",
                  "posterior need not be proper (it usually is not with ",
                  "few groups). The chains still run and Rhat cannot see ",
                  "it. Drop prior = to get the brms default priors ",
                  "instead", call. = FALSE)
    }
    # a MAP fit's penalty is taped INTO fit$obj, so "flat" is a claim
    # about that tape too: the sampled density is the bare likelihood,
    # and the fit's own prior goes with the defaults. Said out loud,
    # because the alternative reading (drop the defaults, keep the
    # fit's penalty) is equally plausible from the argument name
    if (defaults && !is.null(base)) {
      frm_message("prior = \"flat\" drops the prior this fit was made ",
                  "with (", prior_target_list(base), ") as well as the ",
                  "defaults: the sampled density is the bare likelihood, ",
                  "not the penalized objective frm() maximized")
    }
    return(list(effective = NULL, ri = NULL, flat = TRUE))
  }
  st <- prior_stack(base, prior)
  user_pl <- st$pl
  defs <- if (defaults) default_priors_for(fit)
  if (!is.null(defs) && !is.null(user_pl)) {
    user_pl <- keep_default_bounds(user_pl, defs)
    defs <- drop_superseded(defs, user_pl)
  }
  eff <- if (is.null(defs)) user_pl else if (is.null(user_pl)) {
    defs
  } else {
    defs + user_pl
  }
  ri <- if (!is.null(eff)) resolve_prior_input(fit, eff)
  # announced even when there are no defaults at all: a family whose
  # slots frmtmb cannot prior must say so rather than look flat by
  # accident
  if (defaults) {
    announce_default_priors(defs, c(default_prior_notes(fit, user_pl),
                                    unreached_prior_notes(fit, user_pl,
                                                          ri)))
  }
  if (length(st$legacy)) {
    # the legacy spelling takes over exactly the internal parameters it
    # names and leaves the rest of the stack in place
    for (lg in st$legacy) {
      ri <- merge_prior_inputs(ri, resolve_prior_input(fit, lg))
    }
    eff <- eff %||% structure(list(), class = "frmtmb_priorlist")
    attr(eff, "overrides") <- Reduce(utils::modifyList, st$legacy)
  }
  # brms's dirichlet(1) on a mixture's weights is a density on the whole
  # simplex, which no set_prior() class names here (class "theta" is
  # frmtmb's covariance vector), so it rides on the effective prior as
  # a flag that the retape reads (mix_theta_nlp())
  mix <- defaults && length(mix_theta_terms(fit))
  if (mix) {
    eff <- eff %||% structure(list(), class = "frmtmb_priorlist")
    attr(eff, "mix_theta") <- TRUE
  }
  list(effective = eff, ri = ri, mix_theta = mix)
}

#' Sample a model with NUTS
#'
#' Runs [tmbstan::tmbstan()] on the model's objective and returns the
#' draws with frmtmb parameter names. Given a [frmtmb::frm()] fit it samples the
#' fitted objective, initialized at the ML estimates, which shortens
#' warmup considerably. Given a formula and `data` it assembles the same
#' objective without optimizing anything first (see Sampling from a
#' formula).
#'
#' **Both routes sample a posterior under the same default priors** (see
#' Default priors), brms 2.23's own weakly-informative ones. What the
#' fit adds is a starting point, not a different density: the ML mode
#' anchors the chains and shortens warmup, and the model is the one the
#' fit already carries.
#'
#' The two differ only in what else is in the stack. A fit made with
#' `frm(prior = )` is a MAP fit, and it carries its own prior into
#' sampling: the priors that apply are then this call's `prior =` first,
#' the fit's own next, and the defaults under both, most explicit
#' winning per slot it addresses. `prior = "flat"` opts out of all of
#' it and samples the bare likelihood, on either route.
#'
#' [check_laplace()] is the one caller that stays unpenalized by
#' default. It measures the Laplace and Wald approximations against the
#' shape of the very objective the fit maximized, and a default prior
#' would change the thing being measured, so it asks for the flat
#' density explicitly rather than inheriting whatever this function
#' defaults to.
#'
#' @section Sampling from a formula:
#' `frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian())`
#' parses, assembles and tapes exactly as [frmtmb::frm()] does, stops before the
#' optimizer, and hands the objective to Stan. Every pre-optimizer
#' refusal still applies (REML, quadrature, the mixture and `frmtmb.latent::hmm()`
#' guards). There is no mode, so the default `init` is `"random"`:
#' Stan's own overdispersed initialization on the unconstrained scale,
#' inside any bounds the prior carries.
#'
#' The returned object supports the whole draws surface -
#' [summary()], [frmtmb::fixef()], [frmtmb::VarCorr()], [frmtmb::ranef()], [frmtmb::hypothesis()],
#' [posterior_epred()], [posterior_predict()], [posterior_linpred()],
#' [frmtmb::pp_check()], [frmtmb::as_draws()] - because those read the model frame and
#' one draw at a time. The methods that report a maximum-likelihood
#' quantity refuse instead of inventing one: [check_laplace()] (it
#' compares NUTS against a mode that does not exist here), and on the
#' embedded object reachable as `x$fit`, [summary()], [vcov()],
#' [confint()], [logLik()], [frmtmb::fixef()], [frmtmb::ranef()], [frmtmb::VarCorr()],
#' [predict()], [fitted()], [residuals()] and [simulate()].
#'
#' @section Reparameterization:
#' A random-effect block has a FUNNEL in its centered joint posterior:
#' the width of the prior on `b` is a standard deviation that is being
#' sampled at the same time, so the region NUTS must explore narrows as
#' that standard deviation shrinks, and one step size cannot fit both
#' ends of it.
#'
#' `reparameterize = TRUE` (the default) samples `z ~ N(0, I)` instead
#' and computes `b = L(theta) z` on the tape, with `L` the block's own
#' Cholesky factor, which is brms's construction. Each draw is mapped back
#' through ITS OWN `theta`, so the random-effect columns of the draws
#' matrix
#' hold the same quantity in the same order under the same names as
#' `reparameterize = FALSE` gives, and every method downstream
#' ([posterior_epred()], [frmtmb::ranef()], [frmtmb::log_lik()],
#' [frmtmb::loo()], [frmtmb::conditional_effects()],
#' [frmtmb::hypothesis()]) reads them without knowing
#' which route produced them. Only the `stanfit` inside the object
#' carries `z`.
#'
#' Priors are unaffected. They are declared on the OUTER parameters
#' (`beta`, `theta`, `betad`), and `theta` means the same thing on both
#' routes; nothing prior-able is reparameterized. What priors DO decide
#' is which blocks are eligible, below.
#'
#' *Which blocks, and why not all of them.* A block is non-centered when
#' two things hold, and they are the same thing twice: every parameter it
#' has is a standard deviation or a correlation with a Cholesky factor
#' registered for its structure, and every one of those parameters
#' carries a PRIOR. Both are about not handing the sampler a direction
#' it can run away in.
#' Non-centering gives NUTS the run of `theta`'s whole range, and the
#' centered funnel is what was keeping a chain out of the parts of that
#' range where the posterior is flat or improper. Removing the funnel
#' without closing those off first trades a slow chain for a wrong one.
#'
#' In practice that means a default-prior run, on EITHER route, whose
#' priors cover every standard deviation and every correlation,
#' non-centers
#' `(1 | g)` and any block written one term at a time, `diag()` and
#' `homdiag()` blocks, [mgcv::s()] smooths, `equalto()`, `gr(cov = )`,
#' and the CORRELATED blocks `(x | g)`, `cs()`, `ar1()` and `hetar1()`.
#' A `k =` Hilbert-space `gp()` block non-centers too, since its
#' factor is diagonal and the defaults prior its sd and its length
#' scales (brms's classes `"sdgp"` and `"lscale"`), and so does an exact
#' `gp()` field, as brms's `sdgp * L zgp`: its factor is the Cholesky
#' factor of the kernel over the field's distinct positions, taken on the
#' tape at a cost of `n^3 / 3` per gradient for `n` positions.
#' `rr()` is already non-centered by
#' construction, since its own coefficients are the standard normal
#' factors. What stays centered is a run whose variance parameters are
#' flat: `prior = "flat"`, and [check_laplace()], which asks for that
#' density on purpose.
#'
#' The call `message()`s every block it left centered, with the reason.
#' A model made only of those samples centered throughout and says so
#' rather than failing. The reasons:
#'
#' - A FLAT PRIOR on the block's standard deviation. Send `sd` down with
#'   `z` where it is: `b = sd z` goes to zero, the likelihood settles on
#'   the model without that random effect, and the density stops
#'   changing: a flat tail with nothing to stop a chain in it. Measured
#'   on a six-group random-intercept model, one chain of 2000, three
#'   seeds: under flat priors the non-centered chain walks `theta` to
#'   -1e15 at a bulk-ESS of 1; under the default `student_t(3, 0, s)`,
#'   which makes that tail integrable, 174 to 284 against 3 to 48
#'   centered.
#' - A FLAT PRIOR on the block's CORRELATION, which is the same
#'   argument: flat on frmtmb's unbounded correlation parameter is
#'   `(1 - rho^2)^-3/2`, improper, with all its mass at `|rho| = 1`.
#'   Before 0.39 that kept every correlated block centered, because no
#'   proper correlation prior existed to close the tail; `lkj(eta)`
#'   now does, and it is a default, so
#'   `(Days | Subject)`, `cs()`, `ar1()` and `hetar1()` non-center
#'   like any other block.
#' - A Student-t latent (`gr(dist = "student")`): a scale mixture, not a
#'   linear factor.
#' - `car()`, `spde()` and `gr(prec = )`: sparse precisions whose factor
#'   is dense.
#' - The spatial covariances (`ou`, `exp`, `gau`,
#'   `mat`): a full factorization per gradient evaluation.
#' - `toep()`: not positive definite everywhere in its parameterization,
#'   so `b = L z` is not a bijection there.
#'
#' *What it is worth.* One chain of 2000 iterations, three seeds
#' (dev/benchmarks.md). Where each group's own data say little (the
#' regime the funnel lives in) it is decisive: 80 groups of 2 binary
#' observations run at a min-ESS of 236 against 5 centered, 55 times the
#' effective draws per second. Where the groups are informative it is a
#' wash: the `epilepsy` GLMM and an uncorrelated `sleepstudy` are within
#' noise of the centered chain either way. It is not a blanket speed-up,
#' and the blocks it declines to touch are the ones where it would have
#' done harm.
#'
#' `laplace = TRUE` integrates the random effects out, which removes the
#' funnel by itself, so that route ignores `reparameterize` entirely.
#'
#' The ML fit is untouched either way. The Laplace approximation
#' integrates `b` out, and that integral is invariant under a linear
#' change of the integrated variable, so `frm()` has nothing to change
#' and no fitted object ever carries a non-centered tape.
#'
#' @section Laplace draws:
#' `laplace = TRUE` passes to [tmbstan::tmbstan()], which then samples
#' the outer parameters alone (the coefficients, the distributional
#' parameters and the covariance parameters) and integrates the rest
#' out with the Laplace approximation at every step: the group-level
#' coefficients, the coefficients of a smooth, and the values of an
#' `mi()` term. The draws hold no columns for those. brms has no such
#' route, because its sampler always visits the group-level
#' coefficients.
#'
#' *What the draws answer.* Every summary of the outer parameters:
#' `summary()`, `fixef()`, `VarCorr()`, `posterior_summary()`, the
#' `as_draws_*()` converters, `hypothesis()` on the coefficients and on
#' `sd_` and `cor_` quantities, and `mcmc_plot()`. Also every predictive
#' quantity that reads no integrated value: `posterior_epred()`,
#' `posterior_linpred()` and `posterior_predict()` with
#' `re_formula = NA` on a model whose only integrated values are
#' group-level coefficients, a distributional parameter with no
#' group-level term, and the default curves of `conditional_effects()`.
#' These are computed from the draws as they are and are exact.
#'
#' *What they refuse, by name.* A quantity that reads an integrated
#' value: a prediction at the observed group levels (the default
#' `re_formula`, with or without `newdata`), any prediction on a model
#' with a smooth or an `mi()` term that reads them, `ranef()`,
#' `coef()`, `log_lik()` and everything built on it (`loo()`, `waic()`,
#' `psis()`, the `loo_*()` functions, `pp_check()`'s `loo_*` types), and
#' `bayes_R2()`. The
#' refusal names the function you called and says what to do instead;
#' it suggests `re_formula = NA` only where that function takes it and
#' it computes on the model. Whether a call reads an integrated value is
#' found out by evaluating it at one draw with those values set to `NA`
#' and to `0`: a result that changes read them. A draw whose result has
#' non-finite values where the first draw's has none is refused the same
#' way, so a value that vanishes at the first draw only cannot slip
#' through as NaN.
#'
#' *Why they refuse rather than fill the values in.* `laplace = TRUE`
#' samples an approximate posterior, and its use is to check that
#' approximation against draws that sample the group effects. For
#' predictions, and for anything else that conditions on the group
#' effects, sample without `laplace = TRUE`. Filling the values in would
#' also be wrong or random: a prediction at each draw's conditional
#' modes leaves out the uncertainty of the group effects given the
#' outer parameters, and understated the spread of the in-sample
#' `posterior_epred()` by up to a factor of four on a gaussian model
#' with 12 groups of 10 rows; drawing them from their Laplace
#' conditional matched the full draws there, but it would make
#' `posterior_epred()` depend on the random seed and disagree with
#' `log_lik()` and `ranef()` about what a draw contains.
#'
#' *Two edge cases.* On a model with nothing to integrate out,
#' `laplace = TRUE` is the model itself, so the call says so and samples
#' it on the full route, with draws identical to a call without the
#' argument. On a REML fit whose objective is sampled as it stands
#' (`prior = "flat"`), the objective integrates the coefficients out
#' too, and that combination is refused: the draws would hold no
#' coefficient.
#'
#' @section Default priors:
#' Both routes default to brms 2.23's own weakly-informative
#' priors, read off `brms::default_prior()` on matched models, family by
#' family. Write `y*` for the response transformed by a link and
#' `s = max(2.5, round(mad(y*), 1))`:
#'
#' | class | default | scale |
#' |---|---|---|
#' | `b` (slopes) | flat | - |
#' | `Intercept` | `student_t(3, round(median(y*), 1) - m, s)` | link |
#' | `Intercept` (ordinal thresholds) | `student_t(3, 0, 2.5)` | thresholds |
#' | `sd` | `student_t(3, 0, s)` | natural sd, log-Jacobian applied |
#' | `cor`, `rescor`, `lncor`, `cortime` | `lkj(1)` | correlation matrix |
#' | `ar`, `ma` (`cov = TRUE`), `cosy`, `lagsar` | flat, see below | natural |
#' | `simo`, `sbhaz` | `dirichlet(1)` | the simplex |
#' | `theta` (mixture weights) | `dirichlet(1)` | the simplex |
#' | `sigma` | `student_t(3, 0, s)` | natural |
#' | `shape`, `phi` | `gamma(0.01, 0.01)` | natural |
#' | `shape` of the negative binomials | `inv_gamma(0.4, 0.3)` | natural |
#' | `sigma` of `brmsfamily("negbinomial2")`, `1 / shape` | `gamma(0.4, 0.3)` | natural |
#' | `nu` | `gamma(2, 0.1)` | natural |
#' | `kappa` | `gamma(2, 0.01)`; `gamma(0.01, 0.01)` for `xbeta()` | natural |
#' | `beta` (exgaussian) | `gamma(1, 0.1)` | natural |
#' | `zi`, `hu`, `zoi`, `coi`, `bias`, `quantile` | `beta(1, 1)` | natural |
#' | `bs` | `gamma(1, 1)` | natural |
#' | `ndt` | `uniform(0, min(y))` | natural |
#' | `alpha` | `normal(0, 4)` | natural |
#' | `xi` | `normal(0, 2.5)` on `tmp_xi` | see below |
#'
#' The rows marked natural are for a parameter without a formula: the
#' density is on the parameter itself, with the Jacobian of its link,
#' which is what brms's bounds do on its constrained parameters. A
#' parameter with a formula, `shape ~ x` or `shape ~ 1`, has its default
#' on its intercept instead, on the link scale, as in brms:
#' `student_t(3, 0, s)` for `sigma` (with `y*` on `sigma`'s link),
#' `student_t(3, 0, 2.5)` for `shape` and `phi`, `normal(2.7, 0.8)` for
#' `nu`, `normal(5.0, 0.8)` for `kappa`, `normal(1.7, 1.3)` for `beta`,
#' `logistic(0, 1)` for `zi`, `hu`, `zoi`, `coi`, `bias`, `quantile` and
#' a mixture weight, `normal(-0.6, 1.3)` for `bs`, `normal(0, 4)` for
#' `xi` and `alpha`, and `normal(0, 1)` for an ordinal `disc`
#' (`lognormal(0, 1)` under its identity link). `ndt` with a formula
#' stays flat, as in brms. A family brms does not have takes the
#' default of the dpar's name, which is what brms gives a custom family.
#'
#' For `sigma` without a formula `y*` is the response itself (logged for
#' a log-scale family), whatever the `mu` link: brms reads `sigma`'s
#' own link there, which is the identity. `com_poisson()`'s dispersion,
#' frmtmb's `nu`, takes brms's `shape` default, because it is the same
#' quantity under brms's name. A mixture component's dpars take the
#' class's default (`sigma2` is class `sigma`) and no family's own: brms
#' asks the mixture, which states none.
#'
#' The support of each natural default is brms's bound on the class:
#' frmtmb's link maps onto it (`log` onto `(0, Inf)`, `logm1` onto
#' `(1, Inf)` for `nu`, `logit` onto `(0, 1)`). Two classes differ, and
#' the default is still the faithful one. `ndt` of
#' [frmtmb::shifted_lognormal()] is on a log link, and its upper bound,
#' brms's `min_Y`, is the density's: the uniform is `NaN` past it, where
#' the likelihood has no support either. `ndt` of frmtmb.eam's
#' `wiener()` is on a link bounded by the fastest response, so the
#' uniform is flat over exactly brms's interval. A scalar `xi` of
#' [frmtmb::gen_extreme_value()] is sampled through brms's `tmp_xi` (see
#' the message), and is above -1, where brms's range goes on.
#'
#' A residual-correlation coefficient that frmtmb holds through a
#' transform gets brms's `(flat)` on the coefficient itself, written as
#' a uniform density over the whole interval the transform maps onto and
#' carried inward with its Jacobian: `uniform(-1, 1)` on a first-order
#' `cov = TRUE` `ar` or `ma`, a uniform on frmtmb.spatial's `sar()`
#' interval for `lagsar` and `errorsar`, and for `cosy` a uniform on its
#' positive-definite window with `lb = 0`, which is brms's `[0, 1]`. A
#' `cov = TRUE` term of order two or more, which brms does not have,
#' takes a uniform on a box that holds its stationary (or invertible)
#' region, which is flat on that region. A prior of your own on one of
#' these classes keeps the default's bound on each side it leaves
#' unset, and a bounds-only one keeps the flat density, as brms keeps
#' both. The `cov = FALSE` coefficients are unconstrained and are the
#' sampled parameters, so brms's flat needs no density there. A bounded
#' one of these parameters is sampled through one logistic onto its
#' interval, as Stan's bound is in brms.
#'
#' Under a first-order `cov = TRUE` `ar()`, `ma()` or `arma()`, the
#' `sigma` that frmtmb reports (in `summary()`, `sigma()`, `VarCorr()`
#' and these draws) is the MARGINAL sd of a residual, where brms's
#' `sigma` is the innovation sd, smaller by the factor `sqrt(1 - ar^2)`,
#' `1 / sqrt(1 + ma^2)` or `sqrt((1 - ar^2) / (1 + ma^2 + 2 ar ma))`; on
#' an AR(1) of 0.9 that is 1.16 here against 0.43 in brms, for the same
#' model. A prior on class `"sigma"`, the default `student_t(3, 0, s)` or
#' one of your own, is read as brms reads it, on the innovation sd, with
#' the Jacobian of the factor, here and in [frmtmb::frm()] alike, so a
#' `default_prior(route = "sample")` table passed back as `prior =` is
#' the default density. A `lb` above 0 or a `ub` on `sigma` there is
#' refused: a box on the innovation sd is not a box on the marginal sd,
#' which the sampler holds. Under `sigma ~ x` with such a term, a prior on
#' the intercept of `sigma` (`class = "Intercept", dpar = "sigma"`),
#' the default included, is still read at frmtmb's log marginal sd, where
#' brms reads it at the log innovation sd; the slopes are the same in
#' both, as [frmtmb::set_prior()] says.
#'
#' Every default with a bound (`nu` above 1 on an identity link, `sigma`
#' above 0 on one, `ndt` below the smallest response, `cosy` above 0)
#' keeps that bound under a prior of your own on the same slot, on each
#' side the prior leaves unset, as brms keeps the bound it declares.
#'
#' A mixture's weights without a formula are a simplex, and
#' `dirichlet(1)` on it is the flat density there, with the Jacobian of
#' the log ratios the sampler works on. No [frmtmb::set_prior()] class
#' names them (`theta` is frmtmb's covariance vector), so the default
#' is disclosed in the message and `prior = "flat"` drops it with the
#' other defaults.
#'
#' A `mo()` simplex is sampled in softmax coordinates under
#' `dirichlet(1)`, brms's default and the flat density on the simplex,
#' which `prior = "flat"` therefore keeps; the fit's own coordinates
#' (see "Monotonic effects" in [frmtmb::frm()]) have folds at the faces
#' that wall a density off into separate sheets. The draws carry brms's
#' names, `simo_<term>1[k]`, one column per weight. A `simo` row in
#' `prior =`, `set_prior("dirichlet(c(5, 1, 1))", class = "simo",
#' coef = "mox1")`, replaces `dirichlet(1)` on that simplex, as in brms,
#' with the softmax's Jacobian applied.
#'
#' The `disc` default matters more than most: the likelihood cannot
#' tell an intercept in `disc` apart from the scale of the thresholds,
#' so this prior is what places it.
#'
#' `m` is the mean of the predictor's `offset()`, 0 without one. brms
#' takes it off because an offset moves the intercept by that much, and
#' only under the links below that transform the response; under any
#' other link the location stays 0.
#'
#' The `sd` default is written once per response, distributional
#' parameter and nonlinear parameter that has a random-effect standard
#' deviation, as brms writes one row for each, with `y*` on the `mu`
#' link. A class `"sd"` specification reaches only its own response,
#' distributional parameter and nonlinear parameter (see
#' [frmtmb::set_prior()]), so a `set_prior(class = "sd")` of your own
#' replaces the location parameter's default and leaves, for example, a
#' `phi` block's default in place, as in brms.
#'
#' A multivariate model reads `y*` and `s` off each response for that
#' response's rows, which are written with `resp`, as brms's are.
#'
#' The link is transformed only for `identity`, `log`, `inverse`,
#' `sqrt` and `1/mu^2` - brms's own list - with a log-scale family
#' ([frmtmb::lognormal()] and its relatives, whose `mu` link is spelled
#' `identity` but whose response is logged inside the density) counted
#' as `log`. Under any other link (`logit`, `cloglog`, ...), for a
#' response that was a factor, and for a family whose parameters run
#' over categories (ordinal, categorical, multinomial), the location is
#' 0 and the scale 2.5 - which is what brms reports for a bernoulli
#' model. Under a `log`, `inverse` or `1/mu^2` link the ZEROS of the
#' response are replaced by `0.1` before the transform, element by
#' element, exactly as brms does; the non-zero values are left alone,
#' so a count response with median 1 keeps `log(1) = 0` as its
#' location.
#'
#' The call `message()`s the priors it chose, one compact line per
#' class, AND one line per slot it deliberately left flat, so a model
#' that gets few defaults - or none - says so rather than looking flat
#' by accident. Wrap it in `suppressMessages()` to silence that.
#' `prior_summary()` on the returned draws reproduces the chosen priors
#' exactly, except the mixture weights' Dirichlet, which no class names.
#'
#' *What is deliberately NOT matched.* Each of these is named in the
#' message whenever the model has one.
#' - A `toep()` correlation stays flat. Its banded parameterization is
#'   not positive definite everywhere, so there is no correlation matrix
#'   over its whole parameter space for an LKJ density to be about.
#'   Every other correlated structure gets `lkj(1)`, brms's own default
#'   (see [frmtmb::set_prior()] for what the density is on frmtmb's
#'   parameters).
#' - Sum-to-zero ordinal thresholds stay flat: brms puts its density on
#'   thresholds it declares before centering them, a vector frmtmb does
#'   not have.
#' - A Student-t group term's degrees of freedom have no prior: frmtmb
#'   fixes them (`gr(dist_nu = )`), where brms estimates them under
#'   `gamma(2, 0.1)`.
#' - NONLINEAR parameters stay flat, which is what brms does: their
#'   coefficients are class `b`, and the response's median and mad say
#'   nothing about a rate or a shape sitting inside a nonlinear body.
#'   Write them with `set_prior(nlpar = )`, as the brms nonlinear
#'   vignette does.
#'
#' *Overriding and opting out.* A `set_prior()` specification takes over
#' the parameters it reaches. The defaults stay on every other
#' parameter, as in brms. `set_prior(class = "sd", group = "g")`
#' replaces the default on the standard deviations of `g` and keeps it
#' on the other groups. Adding `coef = "x"` narrows the prior to one
#' standard deviation. A default written per coefficient, such as the
#' `lscale` of each `gp()` term, is more specific than a class-wide
#' specification and stays, as in brms, and so does a more specific
#' specification of your own. The default-prior message names a
#' specification that applies to no parameter for either reason.
#' Write the coef to replace such a default. A named list of prior
#' objects takes
#' over exactly the internal parameters it names. `prior = "flat"`
#' turns the defaults off entirely and samples the likelihood, which
#' warns when the model has variance components: their flat-prior
#' posteriors need not be proper, and neither the chains nor Rhat can
#' see that. On a MAP fit `"flat"` reaches further than the defaults:
#' that fit's own prior is taped into its objective, and the bare
#' likelihood is what `"flat"` names, so the penalty is rebuilt out and
#' the call says which prior it dropped.
#'
#' *Three sources, most explicit first.* On a fitted model the stack is
#' this call's `prior =`, then the fit's own prior (`frm(prior = )`
#' makes a MAP fit, and that prior is part of the model the fit is),
#' then the defaults. Each addressed slot is settled by the most
#' explicit source that names it, one slot at a time: a call-level
#' `set_prior(class = "sd")` replaces the fit's `sd` prior and the `sd`
#' default while the `Intercept` default stays. `prior_summary()` on
#' the returned draws prints what the stack came to.
#'
#' @param fit A `frmtmb_fit`, or a `bf()`/formula to assemble and sample
#'   directly (then `data` is required), or a `frmtmb::frm_multiple()`
#'   result (see Multiply imputed data).
#' @param data A data frame of model data, when `fit` is a formula.
#' @param family Family, when `fit` is a plain formula that does not
#'   carry one (`frm_sample(bf(y ~ x), data = dd, family = poisson())`;
#'   the `+` spelling `bf(y ~ x) + poisson()` works too).
#' @param data2,start,fit_control,na.action,REML As in [frmtmb::frm()]
#'   (`fit_control` is that function's `control`). They belong to the
#'   formula path, where they assemble the model this function samples,
#'   and on a fitted object they are refused by name rather than
#'   accepted and discarded: the model is already assembled, so nothing
#'   would read them. `fit_control` is not the sampler's control list;
#'   `control` is.
#' @param control The SAMPLER's control list, brms's spelling and
#'   brms's meaning: `control = list(adapt_delta = 0.99,
#'   max_treedepth = 12)` reaches [tmbstan::tmbstan()] and then rstan
#'   unchanged, so a `brm()` call ports across with this argument
#'   included. `NULL`, the default, leaves rstan's own defaults in
#'   place, and so does `list()`. The names are rstan's (`adapt_delta`,
#'   `max_treedepth`, `stepsize`, `metric`, ...) and an unrecognized
#'   one is refused by name, because rstan answers an unknown option by
#'   returning an empty fit rather than by raising. `frmtmb_control()`
#'   fields arriving here are refused by name too and point at
#'   `fit_control`: earlier releases spelled the fit-time options
#'   `control`, and a silent reinterpretation would hand `grad_tol` to
#'   the sampler. So is an ABBREVIATION of the name: only `control`
#'   binds to this argument, and a shorter spelling would otherwise
#'   travel through `...` and partial-match rstan's own `control`,
#'   reaching the sampler unchecked.
#' @param ... Passed to [tmbstan::tmbstan()] (`chains`, `iter`,
#'   `laplace` (see Laplace draws), `cores`, ...). `cores` parallelizes
#'   over chains on every platform. On Windows the chains run on socket
#'   workers, each of which rebuilds the tape from the serialized
#'   objective closure
#'   (tmbstan retapes on the worker; the closures are self-contained),
#'   giving draws identical to a sequential run at the same seed. The
#'   per-worker startup (a new R process, package load, retape) is a
#'   fixed several seconds, so short chains gain nothing; long chains
#'   approach a per-chain speedup. A core count inherited from
#'   `options(mc.cores)`, which rstan reads when `cores` is not given,
#'   behaves the same way. `seed` fixes Stan's stream, the init
#'   jitter (see `init_jitter`) and the prior draws of
#'   `sample_prior = "yes"`; `NA`, brms's default, leaves all three to
#'   the session. `save_warmup` is rstan's and brms's: `TRUE`, the
#'   default, keeps each chain's warmup iterations beside the draws,
#'   where `inc_warmup = TRUE` reads them (see Warmup).
#' @param sample_prior brms's `sample_prior`: `"no"` (the default),
#'   `"yes"` or `"only"`, or `TRUE` for `"yes"` and `FALSE` for `"no"`.
#'   `"yes"` adds draws from the priors to the draws, under brms's names
#'   (`prior_b`, `prior_Intercept`, `prior_sigma`, ...), which
#'   [prior_draws()] returns and [frmtmb::hypothesis()] reads for the
#'   evidence ratio of a point hypothesis. `"only"` samples the prior of
#'   the whole model with the likelihood left out, so that
#'   [posterior_predict()], [pp_check()] and conditional effects give the
#'   prior predictive distribution. See Prior draws.
#' @param prior Priors: a [frmtmb::set_prior()] specification, or a named list
#'   of prior objects (see [frmtmb::prior_normal()]) whose names are parameter
#'   names as in the draws (or whole components: `"beta"`, `"theta"`,
#'   ...), or the string `"flat"`. The objective is re-taped with the
#'   prior terms added; a fitted model itself is unchanged. On both
#'   paths the brms default priors apply to whatever this argument, and
#'   a MAP fit's own prior, leave alone (see Default priors), and
#'   `prior = "flat"` opts out of them entirely. A `brmsprior` object
#'   built by brms's own `prior()` is translated row by row. The
#'   argument takes brms's spelling, `prior`. On a fit that carries its
#'   own prior, a specification here for a slot the fit's prior also
#'   names (the same class, coef, group, resp, dpar and nlpar) REPLACES
#'   the stored one, as brms's `update(prior = )` replaces a row; it is
#'   not refused as a duplicate, and specifications for other slots
#'   stack with the stored ones, the more specific applying. Two
#'   specifications for one slot within this argument are refused. The
#'   `priors` of releases
#'   before 0.43 is gone rather than aliased, and because this
#'   function's `...` would otherwise swallow it, the old name is
#'   refused by name.
#' @param init Initialization. On a fit the default
#'   (`"last.par.best"`) starts chain 1 exactly at the ML mode and every
#'   further chain at the mode plus a normal perturbation of sd
#'   `init_jitter` on the unconstrained scale. The mode anchor keeps
#'   warmup short; the jitter keeps the chains overdispersed enough for
#'   Rhat to retain power against multimodality (the standard objection
#'   to identical mode starts). `"random"` requests Stan's own
#'   overdispersed initialization, and is the default from a formula,
#'   where there is no mode to start at.
#' @param init_jitter Per-chain perturbation sd for the default init;
#'   `0` starts every chain exactly at the mode. With `seed` given (in
#'   `...`), the jitter is drawn from that seed rather than from the
#'   session's random stream, so `seed` reproduces every chain, as
#'   `brm(seed = )` does. Without one it is drawn from the session's
#'   stream, so `set.seed()` decides it. Either way rstan's sampler
#'   still advances the session's stream (it draws each chain's
#'   permutation of its draws), as it does under `brm()`.
#' @param reparameterize Sample the qualifying random-effect blocks in
#'   their non-centered form (`b = L(theta) z`), which removes the
#'   funnel of the centered joint posterior. `TRUE` by default. Not
#'   every block qualifies, and the call says which did not and why;
#'   see Reparameterization, which also explains why nothing downstream
#'   can tell the two routes apart. The mode-anchored `init` still
#'   starts at the ML mode: it is mapped through
#'   `z0 = L(theta_hat)^-1 b_hat`. Seeded draws differ between the two
#'   routes, so `reparameterize = FALSE` is also the way to reproduce a
#'   run made before this default existed.
#' @param .diagnostic Internal, and named with a dot because there is no
#'   reason to set it by hand. `TRUE` turns the default priors off and
#'   silences the disclosure, and leaves the fit's own objective alone:
#'   a MAP fit's penalty is taped into it, so no prior is resolved
#'   again and no tape is rebuilt. [check_laplace()] sets it, because
#'   it measures the Laplace and Wald approximations against the
#'   density the fit maximized and a default prior would change what is
#'   being measured. Write `prior = "flat"` for the user-facing
#'   opt-out, which goes further and rebuilds a MAP fit's objective
#'   without its penalty.
#' @return An object of class `frmtmb_draws`: list with the `stanfit`,
#'   a draws matrix with named columns (`as.matrix()` method), the
#'   originating fit, and, when any block was non-centered, a
#'   `reparam` note saying which. The `stanfit` holds the parameters as
#'   Stan sampled them, so on a non-centered run its random-effect
#'   columns are `z` while the draws matrix holds `b`. The draws matrix
#'   also carries brms's covariance summaries under brms's names,
#'   `sd_<group>__<coef>`, `cor_<group>__<coef1>__<coef2>`, a smooth's
#'   `sds_`, `rescor__` and the residual structures' `ar[1]`, `cosy` and
#'   `cortime__`, which [frmtmb::variables()] lists for the fit. They
#'   are derived from each draw's internal covariance parameters
#'   (`theta_1`, ...), which stay in the matrix because the model reads
#'   them, so `as_draws_df(draws)$sd_g__Intercept` works as it does on a
#'   brms fit. With `sample_prior = "yes"` the matrix ends with the
#'   prior draws. The object also holds what converts the saved warmup
#'   in the `stanfit` to the draws' columns (and the prior draws of its
#'   rows), the `sample_prior` setting and the arguments of the call,
#'   which [update()][update.frmtmb_draws] reads.
#' @section Prior draws:
#'   `sample_prior = "yes"` draws from the priors as brms does. Each
#'   prior statement of brms's Stan program for the model is one column
#'   of draws, independent of the likelihood: a prior on a whole class is
#'   one column (`prior_b`, `prior_sd_g`, `prior_sigma`); where a
#'   coefficient has a prior of its own, each coefficient of that vector
#'   is a column (`prior_b_x`). A flat prior, the default on class `"b"`,
#'   has no column. A bounded parameter is drawn from its prior truncated
#'   to the bounds: `sigma`, `sd` and `sds` at 0, `nu` at 1, a `"b"` prior
#'   at its `lb` and `ub`. A correlation is one off-diagonal element of an
#'   LKJ draw (`prior_cor_g`, `prior_rescor`), and a `mo()` simplex or a
#'   mixture's weights a Dirichlet draw. `prior_Intercept` is drawn from
#'   the prior of the intercept at the centered predictors, which is the
#'   intercept brms puts that prior on, so the population-level
#'   intercept `b_Intercept` has no prior draws, as in brms. An ordinal
#'   model with a threshold prior of its own (`coef = "1"`) has one
#'   column per threshold, `prior_Intercept__1`, ..., as brms does. A
#'   `gp()` term's `sdgp` and `lscale` are not drawn, nor brms's
#'   `meanme` and `sdme` of a `me()` term, which frmtmb has no prior slot
#'   for, nor a GEV's `xi` (brms draws no prior of it either).
#'
#'   The draws are made directly from each density with R's random
#'   number generators, as brms makes them in Stan's generated
#'   quantities: they are independent, and a run of the sampler on the
#'   prior would give correlated draws of the same distribution at the
#'   cost of a second run. The saved warmup iterations get prior draws
#'   too, as in brms.
#'
#'   `sample_prior = "only"` samples every parameter from its prior with
#'   the likelihood left out (core's `prior_only` objective, documented
#'   in [frmtmb::frmtmb-sampling-api]): the group-level coefficients
#'   from their group's distribution, and the latent values of a
#'   `mi(sdy)` response from their measurement density, as brms does.
#'   The draws are draws from the prior, with the parameter names of the
#'   model, and the predictive methods give the prior predictive
#'   distribution. Every parameter must have a proper prior: a flat one
#'   is refused with brms's message, which names the parameter as brms's
#'   Stan program does (`'b'` with brms's default priors). brms refuses
#'   only a parameter with no prior at all and samples a coefficient
#'   left flat beside others that have a prior; this function refuses
#'   that too. A model with missing values of an `mi()` term is refused,
#'   because without the likelihood those values have no density (brms
#'   samples them flat). A model with `me()` terms is refused: brms
#'   samples it once its classes `"meanme"` and `"sdme"` have proper
#'   priors, and frmtmb has no prior slot for those classes yet. Draws
#'   made with `"only"` have no prior draws of their own, as in brms.
#' @section Warmup:
#'   rstan saves each chain's warmup iterations unless `save_warmup =
#'   FALSE` is passed, and so does brms, whose fit keeps them in its
#'   stanfit. This function keeps them there too, and the
#'   `inc_warmup = TRUE` of [as_draws()] and its relatives, [as.mcmc()],
#'   [rhat()] and [neff_ratio()] moves them to the draws' names and
#'   scales and puts them before each chain's draws, as brms does.
#'   Everything else reads the post-warmup draws only: [ndraws()],
#'   `summary()`, [variables()], the predictions and the diagnostics.
#'   [nuts_params()] and [log_posterior()] read the stanfit, and take
#'   bayesplot's `inc_warmup`.
#' @section The pre-flight against the compatibility registry:
#'   Some models cannot be sampled at all, and the reason usually
#'   belongs to the package that supplies the feature rather than to
#'   this one. Such a package declares the pair through
#'   [frmtmb::frmtmb_register_compat()] with the status `"refused"`,
#'   and this function reads those rows before it does anything else:
#'   before the rstan and tmbstan namespaces are loaded, before the
#'   formula route calls [frmtmb::frm()], and before either route
#'   retapes the objective for its priors. [as_tmbstan()] runs the same
#'   check. The refusal repeats the registry's own note and names the
#'   row, so `frm_compat("frm_sample", "<feature>")` gives the same
#'   answer without running anything.
#'
#'   *What it can decide, and what it will not guess at.* Whether a
#'   model uses a registry feature is settled by WHERE a call sits in
#'   the formula. The check decides three ways.
#'
#'   A response FAMILY is matched by name, which is exact, because a
#'   family is named and never called. It does not look inside
#'   [frmtmb::mixture()], whose family name is the whole string
#'   `"mixture(gaussian, gaussian)"`, so a refused `gaussian` row does
#'   not reach a mixture of gaussians.
#'
#'   An ADDITION TERM is matched by POSITION. It is written inside the
#'   bar on the left of the response formula and nowhere else, and the
#'   `bf()` object carries that, so `y | se(v) ~ x` is separated from
#'   `y ~ a * se(x)` with a user's own `se()`, which frmtmb fits.
#'
#'   Everything else the vocabulary writes as a call is matched BY
#'   NAME, on the display name, outside addition-term position, and
#'   only when a formula actually writes that name and it is not also a
#'   base R function. This is the part that is justifiable from names
#'   alone rather than sound: a name test cannot see that a call sits
#'   where the feature cannot occur, so a refusal registered on `mo()`,
#'   `ma()`, `cosy()` or `mm()` would also fire on a user's own
#'   function of that name in a nonlinear body. Nothing is registered
#'   on those today. `mi_pred()`, `gp_pred()` and `cs_pred()` are
#'   spelled in the vocabulary differently from the way a formula
#'   writes them, which the registry itself reports, so they warn
#'   rather than matching a call no formula makes.
#'
#'   A refused row naming anything else, a covariance structure, a
#'   fitting mode, another post-fit method, WARNS and samples, and so
#'   does one whose name is a base R function. Refusing would fire on
#'   correct models; silence would make a registered refusal look like
#'   a passed check.
#'
#'   `frmtmb.ode` registers the one refused row there is: an
#'   `frm_ode()` fit aborts at the first warmup step inside RTMBode,
#'   and the abort leaves rstan's autodiff arena unusable for the rest
#'   of the session, so the refusal is worth more than the attempt.
#' @section Multiply imputed data:
#'   Given a `frmtmb::frm_multiple()` result, `frm_sample()` samples the
#'   fit of every imputation with the same arguments and combines the
#'   draws, as brms's `brm_multiple()` combines the chains of its fits:
#'   the draws of imputation 2 follow those of imputation 1, and so on,
#'   and the `stanfit` holds every chain, so `nchains()` is the number
#'   of imputations times `chains`. A `seed` gives imputation `j` the
#'   seed `seed + j - 1`, so the result is reproducible and no two
#'   imputations share a random stream. The draws object keeps the
#'   first imputation's fit, so a method that reads the data in sample
#'   (`posterior_predict()` without `newdata`, `pp_check()`, `loo()`)
#'   reads the first imputed data set, which is what brms does and warns
#'   about. Chains of different imputations sample different
#'   posteriors, so an R-hat across them measures the imputations as
#'   well as the sampler, as the brms vignette on missing values says;
#'   check each imputation's chains on their own.
#' @section Multimodal posteriors:
#'   For [frmtmb::mixture()] fits the posterior is multimodal by construction
#'   (label switching at minimum). Mode-centered inits, jittered or
#'   not, leave every chain in one symmetry branch, so Rhat cannot
#'   flag the others; use `init = "random"` there and inspect chains
#'   individually.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#' set.seed(9)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)
#' summary(ds)
#' fixef(ds)
#' hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2) = 0",
#'            class = NULL)
#'
#' # the same model sampled straight from the formula, with no ML fit.
#' # It reports the brms default priors it chose, and prior_summary()
#' # gives them back.
#' ds2 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
#'                   chains = 1, iter = 500, refresh = 0)
#' prior_summary(ds2)
#' fixef(ds2)
#'
#' # a set_prior() specification takes over the parameters it reaches
#' # and leaves the defaults on the rest
#' ds3 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
#'                   chains = 1, iter = 500, refresh = 0,
#'                   prior = set_prior("exponential(1)", class = "sd"))
#' prior_summary(ds3)
#' }
#' }
#' @export
frm_sample <- function(fit, data = NULL, family = NULL, ...,
                       prior = NULL, sample_prior = "no", init = NULL,
                       init_jitter = 0.25, reparameterize = TRUE,
                       data2 = list(), start = NULL,
                       control = NULL,
                       fit_control = frmtmb_control(),
                       na.action = stats::na.omit, REML = FALSE,
                       .diagnostic = FALSE) {
  if (inherits(fit, "frmtmb_multiple")) {
    return(sample_multiple(fit, match.call(), parent.frame()))
  }
  # `...` goes to tmbstan, which would take the retired spelling as an
  # unknown sampler option and run the model with no priors at all
  refuse_retired_priors(list(...), "frm_sample()")
  sample_prior <- check_sample_prior(sample_prior)
  if (isTRUE(.diagnostic) && !identical(sample_prior, "no")) {
    frm_stop("check_laplace() compares the sampler with the density the ",
             "fit maximized, which has no prior draws and keeps the ",
             "likelihood; it cannot take sample_prior = \"", sample_prior,
             "\"", call. = FALSE)
  }
  # what update() needs to run this call again with one argument changed
  sample_args <- list(prior = prior, sample_prior = sample_prior,
                      init = init, init_jitter = init_jitter,
                      reparameterize = reparameterize, control = control)
  sample_args <- c(sample_args[!vapply(sample_args, is.null, NA)],
                   list(...))
  # `control` binds only on its exact spelling here, but an
  # abbreviation is not therefore harmless: it lands in `...` and
  # partial-matches rstan's own `control` further down
  refuse_partial_control(list(...), "frm_sample()")
  check_stan_control(control, "frm_sample()")
  from_formula <- !inherits(fit, "frmtmb_fit")
  if (!from_formula) {
    if (!is.null(data) || !is.null(family)) {
      frm_stop("frm_sample(data =, family =) belongs to the formula ",
               "interface; the model of a fitted object is already fixed. ",
               "Drop them, or pass the formula instead of the fit",
               call. = FALSE)
    }
    # The same rule, applied to the rest of the assembly arguments,
    # which used to be accepted and discarded. `fit_control` is the one
    # that matters most: it is the name the `control` refusal points a
    # confused user at, so accepting a sampler list under it and
    # sampling with the defaults would recreate, under a new name, the
    # silence this release removes. missing() rather than a default
    # comparison, so a caller who passes the default explicitly is
    # still told the argument does nothing here.
    inert <- c(if (!missing(fit_control)) "fit_control",
               if (!missing(start)) "start",
               if (!missing(data2)) "data2",
               if (!missing(na.action)) "na.action",
               if (!missing(REML)) "REML")
    if (length(inert)) {
      frm_stop("frm_sample(", paste0(inert, " =", collapse = ", "),
               ") belongs to the formula interface: those arguments ",
               "assemble a model, and a fitted object is already ",
               "assembled, so they were read by nothing. Pass the formula ",
               "instead of the fit, or refit with frm(). The SAMPLER's ",
               "options are `control`, control = list(adapt_delta = 0.99)",
               call. = FALSE)
    }
  }
  # Before the namespace loads, before frm() on the formula route and
  # before any retape on the fit route, so a refused model costs the
  # registry read and nothing else. The pre-flight needs neither rstan
  # nor tmbstan, and a model that cannot be sampled at all is a more
  # useful thing to be told than which package is missing.
  sample_preflight(fit, family)
  prior <- as_priorlist(prior)
  # the user's own specification, before it is stacked on the fit's
  # prior and the defaults, which may legitimately share a slot
  check_prior_slots(prior)
  if (!requireNamespace("tmbstan", quietly = TRUE) ||
      !requireNamespace("rstan", quietly = TRUE)) {
    frm_stop("frm_sample() needs the 'tmbstan' and 'rstan' packages",
             call. = FALSE)
  }
  if (!is.logical(reparameterize) || length(reparameterize) != 1L ||
      is.na(reparameterize)) {
    frm_stop("reparameterize = must be TRUE or FALSE: it selects the ",
             "non-centered sampling parameterization of the random-effect ",
             "blocks that have one", call. = FALSE)
  }
  if (from_formula) {
    fit <- sample_assemble(fit, data, family, data2 = data2,
                           start = start, fit_control = fit_control,
                           na.action = na.action, REML = REML)
  }
  # no mode to anchor on when the model was never optimized
  init <- init %||% if (from_formula) "random" else "last.par.best"
  pr_lower <- c()
  pr_upper <- c()
  obj <- fit$obj
  # BOTH routes sample a POSTERIOR, so both default to brms's
  # weakly-informative priors: a literal translation of a prior-free
  # brm() call onto a fitted model otherwise pins the chain at whatever
  # boundary mode maximum likelihood found. Precedence is
  # most-explicit-wins per slot - this call's prior =, then the fit's
  # own (a MAP fit's), then the defaults - and prior = "flat" opts out
  # of all of it.
  #
  # `own` is read before it is overwritten, and `.diagnostic` leaves it
  # out of the stack on purpose: a MAP fit's penalty is taped INTO
  # fit$obj at fit time, so resolving it again here would rebuild a
  # tape that check_laplace() already has. See check_laplace().
  own <- fit[["prior"]]
  rp <- sample_resolve_priors(fit, prior,
                              base = if (!isTRUE(.diagnostic)) own,
                              defaults = !isTRUE(.diagnostic))
  fit$prior <- rp$effective
  ri <- rp$ri
  if (!is.null(ri)) {
    pr_lower <- ri$lower
    pr_upper <- ri$upper
  }
  # brms's prior statements of this model, which "yes" draws from and
  # "only" requires to be proper
  plan <- if (!identical(sample_prior, "no")) {
    sample_prior_plan(fit, rp, sample_prior)
  }
  prior_only <- identical(sample_prior, "only")
  if (prior_only) {
    prior_only_check(plan)
    # a missing value of an mi() term has no density without the
    # likelihood; a mi(sdy) response's latent values keep their
    # measurement density, as brms keeps it
    mim <- fit$frame[["mi_map"]] %||% list()
    nmiss <- vapply(mim, function(mm) {
      length(setdiff(mm[["rows"]], mm[["obs"]] %||% integer(0)))
    }, 1L)
    if (any(nmiss > 0L)) {
      frm_stop("frm_sample(sample_prior = \"only\") cannot sample this ",
               "model: ", paste0(nmiss[nmiss > 0L], " missing value",
                                ifelse(nmiss[nmiss > 0L] == 1L, "", "s"),
                                " of ", names(mim)[nmiss > 0L],
                                collapse = " and "),
               " have no density without the likelihood, so their prior ",
               "is flat and improper. brms samples them flat; ",
               "frm_sample() refuses", call. = FALSE)
    }
    # brms needs proper priors on meanme and sdme here, and refuses with
    # its first me() term's meanme otherwise; frmtmb has no slot for
    # either class
    if (!is.null(fit$frame[["me"]])) {
      prior_only_stop("meanme_1", paste0(
        "frmtmb has no prior slot yet for brms's classes \"meanme\" and ",
        "\"sdme\", the mean and SD of a me() term's latent values, so a ",
        "model with me() terms cannot be given the proper priors brms ",
        "samples it from"))
    }
    unc <- prior_only_uncovered(fit, ri$entries, isTRUE(rp$mix_theta))
    if (length(unc)) prior_only_stop(unc[1L])
  }
  # a bounded residual-correlation parameter is sampled through one
  # logistic onto its interval (chart_fn()), not through Stan's bound on
  # the internal parameter, which the chart then replaces. Not under
  # .diagnostic, which may sample fit$obj as it stands
  chart <- if (!isTRUE(.diagnostic) && !is.null(ri)) {
    autocor_bound_chart(fit, ri$lower, ri$upper)
  }
  if (length(chart)) {
    cn <- vapply(chart, `[[`, "", "name")
    pr_lower <- pr_lower[setdiff(names(pr_lower), cn)]
    pr_upper <- pr_upper[setdiff(names(pr_upper), cn)]
  }
  # set_prior()'s lb/ub is the only source of bounds: frm_sample() has
  # no lower/upper of its own, and every outer parameter Stan can
  # constrain has a class that addresses it
  lower <- as.list(pr_lower)
  upper <- as.list(pr_upper)
  laplace <- isTRUE(list(...)$laplace)
  drop_laplace <- laplace && !length(intersect(
    c("b", "miss"), names(fit$frame[["par_template"]])))
  if (drop_laplace) {
    # tmbstan negates obj$env$random, which is NULL here, and died in
    # R's own "invalid argument to unary operator". A Laplace
    # approximation over nothing is the model itself, so the full route
    # samples exactly the density that was asked for
    frm_message("frm_sample(laplace = TRUE): the model has no random ",
                "effects or mi() values to integrate out, so the Laplace ",
                "route and the full route sample the same density. ",
                "Sampling the model as it is")
    laplace <- FALSE
  }
  # frm() skips this check for the unfitted object, because only here
  # are the priors known that the chains will sample under. brms's
  # default on sigma is finite at 0, so with weights summing to n + 1
  # or more over rows that each own an effect, the posterior does not
  # integrate and the chains do not mix (Rhat 2.5 on 2 chains). The fit
  # route was refused by frm() already, or fitted under a prior that
  # bounds it
  if (from_formula && !prior_only) {
    check_weights_unbounded(fit$spec, fit$frame, rp$effective,
                            exact = !laplace, route = "sample")
  }
  # the priors ride on the OUTER parameters, so they compose with the
  # reparameterization rather than interacting with it: theta is theta on
  # either route, and no prior addresses b (or z) at all
  # a scalar GEV xi is sampled through brms's tmp_xi where it can be,
  # and said so where it cannot (gev_xi_blocked_text())
  gev_terms <- gev_xi_chart_terms(fit)
  gev_chart <- any(vapply(gev_terms, function(tm) !length(tm$blocked), NA))
  for (tm in gev_terms) {
    if (length(tm$blocked)) {
      frm_warning("frm_sample(): ", gev_xi_blocked_text(tm, ""),
                  call. = FALSE)
    }
  }
  ncp <- ncp_plan(fit, reparameterize, laplace, ri$entries)
  announce_ncp(ncp)
  # every prior-carrying route rebuilds from build_objective(), the
  # BARE likelihood, and adds the resolved entries once; fit$obj is
  # sampled untouched only when nothing has to be added to it. Under
  # sample_prior = "only" the tape is core's prior_only one, built from
  # a copy of the frame, so the stored fit keeps its likelihood
  tfit <- fit
  if (prior_only) tfit$frame[["prior_only"]] <- TRUE
  if (length(ncp$idx)) {
    obj <- ncp_objective(tfit, ncp$idx, ri$entries %||% list(), chart)
  } else if ((!is.null(ri) && length(ri$entries)) || isTRUE(rp$mix_theta) ||
               length(chart) || prior_only) {
    obj <- prior_augmented_obj(tfit, ri$entries %||% list(), chart)
  } else if (isTRUE(rp$flat) && !isTRUE(.diagnostic) && !is.null(own)) {
    # the opt-out has something to SUBTRACT here: a MAP fit's own prior
    # is part of fit$obj, and "flat" on the user-facing argument means
    # the likelihood without it. Under `.diagnostic` it means the
    # weaker thing - no prior ADDED to the density the fit maximized -
    # which is fit$obj as it stands
    obj <- prior_augmented_obj(fit, list())
  } else if (length(simplex_frame_terms(fit)) || gev_chart) {
    # fit$obj reads a mo() simplex in the maximum-likelihood chart, which
    # has no density a sampler can use; the retape samples the softmax
    # under dirichlet(1). A MAP fit under `.diagnostic` (check_laplace())
    # keeps the prior it was fitted with, which is taped into fit$obj and
    # is resolved again for the retape
    own_entries <- if (isTRUE(.diagnostic) && !is.null(own)) {
      sample_resolve_priors(fit, NULL, base = own,
                            defaults = FALSE)$ri$entries
    }
    obj <- prior_augmented_obj(fit, own_entries %||% list())
  }
  # the sampled simplex coordinates are softmax ones exactly when the
  # objective was retaped (mo_sample_frame())
  mo_soft <- length(simplex_frame_terms(fit)) && !identical(obj, fit$obj)
  if (laplace && "beta" %in% names(obj$env$par)[obj$env$random]) {
    # a REML fit's own objective integrates the coefficients out too.
    # The prior-carrying routes rebuild it without that, but sampled as
    # it stands (prior = "flat", check_laplace()) tmbstan's Laplace route
    # draws theta and sigma alone, and the labels below, which assume
    # the coefficients are there, would name those draws b_Intercept
    frm_stop("frm_sample(laplace = TRUE) on this REML fit would integrate ",
             "the population-level coefficients out along with the random ",
             "effects, and the draws would hold no coefficient at all. ",
             "Sample a maximum-likelihood fit (REML = FALSE) with ",
             "laplace = TRUE, or this fit without it", call. = FALSE)
  }
  bounds <- if (length(lower) || length(upper)) {
    resolve_bounds(fit, unlist(lower), unlist(upper))
  }
  if (identical(init, "last.par.best")) {
    # a singular ML mode (variance at the boundary) is exactly the
    # pathological start the mode-init criticism is about. Only the
    # log-sd components read that way: a CAR mixing proportion or an
    # AR(1) phi at its own boundary is a large theta on a converged,
    # perfectly samplable fit (see log_sd_theta_index()).
    th_all <- fit$estimates[["theta"]] %||% numeric(0)
    sd_i <- log_sd_theta_index(fit)
    ext <- sd_i[abs(th_all[sd_i]) > 8]
    if (length(ext)) {
      frm_warning("The ML mode has an extreme covariance parameter (",
                  paste(names(ext), collapse = ", "),
                  "; likely a boundary/singular fit); mode initialization ",
                  "starts the chains there. The default priors usually ",
                  "pull a chain off that boundary; if this one stays ",
                  "pinned, use init = \"random\" or a tighter prior =",
                  call. = FALSE)
    }
    chart_mode_warning(fit, chart)
    lpb <- obj$env$last.par.best
    rnd <- obj$env$random
    # under laplace tmbstan samples only the outer parameters, so the
    # full-length mode has the wrong length; take the outer slice
    mode <- if (laplace && length(rnd)) lpb[-rnd] else lpb
    if (!is.null(bounds)) {
      # a bound tighter than the ML estimate makes the mode inadmissible
      # as a start; the chains then begin somewhere else, and a
      # posterior pinned against a bound is a modeling statement worth
      # hearing about
      om <- if (laplace || !length(rnd)) mode else mode[-rnd]
      viol <- which(as.numeric(om) < bounds$lower |
                      as.numeric(om) > bounds$upper)
      if (length(viol)) {
        frm_warning("The ML mode violates the requested bound(s) on ",
                    paste(outer_par_names(fit)[viol], collapse = ", "),
                    "; every chain starts at the clamped value instead of ",
                    "the mode. The bounded posterior is not centered on the ",
                    "unconstrained ML estimate", call. = FALSE)
      }
    }
    mb <- mode_aligned_bounds(obj, bounds, laplace, length(mode))
    # the jitter is the one draw R makes for the sampler; with a seed it
    # comes from that seed, so `seed =` fixes every chain as it does in
    # brms, where Stan draws the inits
    init <- with_sample_seed(list(...)$seed, mode_inits(
      mode, list(...)$chains %||% 4, init_jitter, mb$lower, mb$upper))
  }
  args <- list(obj = obj, init = stan_init_arrays(init), ...)
  if (drop_laplace) args$laplace <- NULL
  # tmbstan forwards `control` to rstan::sampling() untouched
  if (!is.null(control)) args$control <- control
  # A development namespace cannot travel to a worker, so this is caught
  # here by name rather than reported by rstan as empty chains.
  dev <- if (.Platform$OS.type == "windows" && stan_cores(args) > 1 &&
               (args$chains %||% 4) > 1) dev_namespaces_of(fit) else character(0)
  if (length(dev)) {
    frm_stop("frm_sample(): cores = ", stan_cores(args), " starts one R ",
             "process per chain on Windows, and each rebuilds the model from ",
             "the packages in the library. This model's family comes from ",
             paste(dev, collapse = " and "), " as loaded by ",
             "pkgload::load_all(), which exists only in this session: a ",
             "worker replaces it with the global environment and fails at the ",
             "family's first internal call. Install the package, or run with ",
             "cores = 1", call. = FALSE)
  }
  # rstan runs parallel chains on PSOCK workers on Windows. The tape's
  # external pointer dies in serialization, but tmbstan ships the cure:
  # its sampling method calls fn() on the worker, which retapes from
  # the objective closure, and frmtmb's generated closures are
  # self-contained (data baked in, namespace references loading frmtmb
  # and RTMB on deserialization), so each worker rebuilds its own tape.
  # Verified seed-for-seed identical to sequential chains. The only
  # cost is startup, hence the note rather than a fallback.
  if (.Platform$OS.type == "windows" && stan_cores(args) > 1 &&
      (args$chains %||% 4) > 1) {
    frm_message("parallel chains on Windows start one R process per core ",
                "and rebuild the tape in each: expect several seconds of ",
                "startup before sampling. Short chains may run faster ",
                "with cores = 1")
  }
  if (!is.null(bounds)) {
    args$lower <- bounds$lower
    args$upper <- bounds$upper
  }
  check_tmbstan_build("frm_sample()")
  sf <- do.call(tmbstan::tmbstan, args)
  check_stan_draws(sf, "frm_sample()")
  # laplace draws skip the inner components entirely; labeling them
  # with the full template order would misattribute theta as a random
  # effect. brms's names where brms has the parameter (b_, r_), the
  # internal ones where it does not (theta_1): see brms_par_labels()
  # how the sampled coordinates become the draws: the same steps for the
  # draws and, when inc_warmup asks, for the warmup (draws_warmup())
  recipe <- list(laplace = laplace, mo_soft = mo_soft, ncp_idx = ncp$idx,
                 chart_lower = if (length(chart)) ri$lower,
                 chart_upper = if (length(chart)) ri$upper)
  gev <- if (gev_chart && !identical(obj, fit$obj)) {
    list(obj = obj, terms = gev_terms)
  }
  m <- sample_finish(rstan::extract(sf, permuted = FALSE), fit, recipe,
                     chart, gev)
  # the warmup rstan saved (save_warmup = TRUE, its default and brms's)
  # stays in the stanfit alone, as brms keeps it, and inc_warmup = TRUE
  # converts it when asked (draws_warmup()). A GEV xi chart reads every
  # draw through the objective's tape, which the draws do not keep, so
  # there the warmup is converted now
  nw <- sf@sim[["warmup2"]][1L] %||% 0
  n_warm <- if (isTRUE(nw > 0)) nw * (sf@sim[["chains"]] %||% 1L) else 0
  warm <- if (!is.null(gev) && isTRUE(nw > 0)) {
    aw <- rstan::extract(sf, permuted = FALSE, inc_warmup = TRUE)
    sample_finish(aw[seq_len(nw), , , drop = FALSE], fit, recipe, chart,
                  gev)
  }
  warm_prior <- NULL
  pinfo <- NULL
  if (identical(sample_prior, "yes")) {
    pd <- with_sample_seed(list(...)$seed, prior_draws_make(
      plan, nrow(m) + n_warm, fit$prior))
    if (!is.null(pd)) {
      # the post-warmup rows first, so that they do not depend on
      # whether the warmup was saved
      m <- cbind(m, pd[seq_len(nrow(m)), , drop = FALSE])
      if (n_warm > 0) {
        warm_prior <- pd[nrow(m) + seq_len(n_warm), , drop = FALSE]
        if (!is.null(warm)) {
          warm <- cbind(warm, warm_prior)
          warm_prior <- NULL
        }
      }
    }
  }
  if (!identical(sample_prior, "no")) {
    pinfo <- list(intercepts = attr(plan, "intercepts"))
  }
  sample_args$own <- list(own)
  sample_args$from_formula <- from_formula
  structure(list(stanfit = sf, draws = m, fit = fit,
                 reparam = if (length(ncp$idx)) {
                   list(blocks = ncp$idx, labels = ncp$labels,
                        centered = ncp$centered)
                 },
                 warmup = warm, warmup_recipe = recipe,
                 warmup_prior = warm_prior, sample_prior = sample_prior,
                 prior_info = pinfo, sample_args = sample_args),
            class = "frmtmb_draws")
}

#' An iterations x chains x parameters array of sampled coordinates as
#' the chain-major draws matrix every method reads: brms's names, the
#' residual-correlation charts, a GEV `xi` (`gev`, the objective and its
#' terms), the `mo()` simplexes and the non-centered blocks moved back,
#' and the unmodeled distributional parameters on their own scale.
#' `recipe` is what `frm_sample()` stores for the warmup; `chart` is
#' rebuilt from it when `NULL` is not given.
#'
#' @noRd
sample_finish <- function(a, fit, recipe, chart = NULL, gev = NULL) {
  laplace <- isTRUE(recipe$laplace)
  # laplace draws skip the inner components entirely; labeling them
  # with the full template order would misattribute theta as a random
  # effect. brms's names where brms has the parameter (b_, r_), the
  # internal ones where it does not (theta_1): see brms_par_labels()
  labels <- brms_par_labels(fit, include_random = !laplace)
  stan_names <- dimnames(a)[[3]]
  n_lab <- min(length(labels), length(stan_names))
  stan_names[seq_len(n_lab)] <- labels[seq_len(n_lab)]
  m <- do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch) a[, ch, ]))
  colnames(m) <- stan_names
  # first, while every column is still in the sampled coordinates the
  # objective reads
  m <- chart_draws_back(m, fit, chart, laplace)
  if (!is.null(gev)) m <- gev_draws_to_xi(m, gev$obj, fit, gev$terms, laplace)
  if (isTRUE(recipe$mo_soft)) m <- mo_draws_to_chart(m, fit, laplace)
  if (length(recipe$ncp_idx)) m <- ncp_backtransform(m, fit, recipe$ncp_idx)
  # brms reports an unmodeled distributional parameter as the parameter
  # itself: the `sigma` column holds sigma, not log sigma
  draws_to_natural(m, fit)
}

#' The saved warmup of draws in the draws' coordinates, or `NULL` when
#' rstan saved none. `frm_sample()` keeps it in the stanfit alone and
#' converts it here, with the recipe of its draws and the prior draws it
#' made for the warmup rows; draws that carry a converted `warmup` (a
#' GEV xi chart, several imputations) return it.
#'
#' @noRd
draws_warmup <- function(x) {
  if (!is.null(x$warmup)) return(x$warmup)
  rc <- x$warmup_recipe
  sf <- x$stanfit
  if (is.null(rc) || !inherits(sf, "stanfit")) return(NULL)
  nw <- sf@sim[["warmup2"]][1L] %||% 0
  if (!isTRUE(nw > 0)) return(NULL)
  chart <- if (!is.null(rc$chart_lower) || !is.null(rc$chart_upper)) {
    autocor_bound_chart(x$fit, rc$chart_lower, rc$chart_upper)
  }
  aw <- rstan::extract(sf, permuted = FALSE, inc_warmup = TRUE)
  w <- sample_finish(aw[seq_len(nw), , , drop = FALSE], x$fit, rc, chart)
  if (!is.null(x$warmup_prior)) w <- cbind(w, x$warmup_prior)
  w
}

#' The prior statements of the model `fit` holds under its effective
#' prior: `validate_prior()`'s table of that prior, read as brms's Stan
#' program writes it (`prior_plan()`). The named-list prior spelling
#' takes over internal parameters that no `set_prior()` row describes,
#' so brms's draws of the rows would not be the prior that was sampled
#' under; it is refused.
#'
#' @noRd
sample_prior_plan <- function(fit, rp, sample_prior) {
  eff <- rp$effective
  if (length(attr(eff, "overrides") %||% list())) {
    frm_stop("frm_sample(sample_prior = \"", sample_prior, "\") reads the ",
             "prior as set_prior() rows, as brms does, and the named-list ",
             "prior spelling takes over internal parameters no row ",
             "describes. Write the prior with set_prior()", call. = FALSE)
  }
  pl <- if (length(unclass(eff %||% list()))) {
    structure(unclass(eff), class = "frmtmb_priorlist")
  }
  tab <- validate_prior(pl, fit, route = "fit")
  # the posterior draws' names, which give brms's spelling of the
  # simplexes, smooths and correlations
  nc <- draws_natural_cols(fit)
  cn <- unique(c(brms_par_labels(fit, include_random = FALSE), nc$extra,
                 nc$derived$names,
                 unlist(lapply(nc$ordinal, `[[`, "names"))))
  prior_plan(fit, tab, cn, mix_theta = isTRUE(rp$mix_theta))
}

#' `frm_sample()` of a `frm_multiple()` result: each imputation's fit
#' sampled by the same call, and the draws combined as brms's
#' `brm_multiple()` combines its chains (`rstan::sflist2stanfit()`).
#' `cl` is the caller's call, re-evaluated once per imputation in
#' `env` with that imputation's fit in place of the pooled object, so
#' an argument the caller did not give stays missing there and the
#' refusals of the fit route still read it so.
#'
#' @noRd
sample_multiple <- function(x, cl, env) {
  fits <- x[["fits"]]
  seed <- if (!is.null(cl[["seed"]])) eval(cl[["seed"]], env)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                           is.na(seed))) {
    frm_stop("frm_sample(): `seed` must be one number, not ",
             deparse1(seed), call. = FALSE)
  }
  ds <- lapply(seq_along(fits), function(j) {
    cj <- cl
    # the function itself, so the call needs no frm_sample on the path
    cj[[1L]] <- frm_sample
    cj[["fit"]] <- fits[[j]]
    # a random stream per imputation, reproducible from the one seed
    if (!is.null(seed)) cj[["seed"]] <- seed + j - 1
    eval(cj, env)
  })
  nm <- colnames(ds[[1L]]$draws)
  for (j in seq_along(ds)[-1L]) {
    if (!identical(colnames(ds[[j]]$draws), nm)) {
      frm_stop("frm_sample() of a frm_multiple() result combines the ",
               "draws of every imputation, and imputation ", j, " has ",
               "different parameters from the first: a factor or group ",
               "level, say, is missing from one of the imputed data ",
               "sets. Sample each fit on its own, x$fits[[j]]",
               call. = FALSE)
    }
  }
  out <- ds[[1L]]
  out$stanfit <- rstan::sflist2stanfit(lapply(ds, `[[`, "stanfit"))
  out$draws <- do.call(rbind, lapply(ds, `[[`, "draws"))
  # the warmup in the stanfit's chain order, each imputation's converted
  # with its own fit, since the combined draws keep the first fit alone;
  # NULL when none was saved
  ws <- lapply(ds, draws_warmup)
  out$warmup_recipe <- NULL
  out$warmup_prior <- NULL
  if (all(vapply(ws, function(w) !is.null(w), NA))) {
    out$warmup <- do.call(rbind, ws)
  } else {
    out$warmup <- NULL
  }
  out$nimp <- length(ds)
  out
}

#' @rdname sample-as_draws
#' @export
as.matrix.frmtmb_draws <- function(x, pars = NA, variable = NULL,
                                   draw = NULL, subset = NULL, ...) {
  # brms:::as.matrix.brmsfit: an unclassed draws_matrix, so the margins
  # are named `draw` and `variable` and the chain count rides along.
  # `x$draws` is the plain matrix, for anyone who wants that
  unclass(draws_accessor_args(x, pars, variable, draw, subset,
                              "as.matrix()",
                              format = posterior::as_draws_matrix, ...))
}

#' @export
print.frmtmb_draws <- function(x, ...) {
  frm_check_dots(...)
  m <- x$draws
  keep <- draws_outer_cols(x)
  cat("<frmtmb_draws> ", nrow(m), " draws x ", ncol(m),
      " parameters\n\n", sep = "")
  tab <- t(vapply(keep, function(nm) {
    c(mean = mean(m[, nm]), sd = stats::sd(m[, nm]),
      `2.5%` = unname(stats::quantile(m[, nm], 0.025)),
      `97.5%` = unname(stats::quantile(m[, nm], 0.975)))
  }, numeric(4)))
  print(signif(tab, 4))
  invisible(x)
}

#' Check the Laplace/Wald approximation against NUTS
#'
#' Samples the fitted objective (see [frm_sample()]) and compares the ML
#' estimates and sdreport standard errors against posterior means and
#' SDs. Close agreement supports the Laplace approximation and Wald
#' intervals; a posterior SD much larger than the Wald SE, or a shifted
#' mean, flags parameters where they are unreliable (typically variance
#' components with few groups).
#'
#' **It samples the density the fit maximized, and no other.** It hands
#' Stan `fit$obj` as it stands. For an ordinary fit that objective is
#' the LIKELIHOOD, so the run is flat; for a MAP fit (`frm(prior = )`)
#' the penalty is taped INTO that objective, so the run carries it,
#' which is right: the mode and the Wald standard errors this function
#' compares against come from the same penalized Hessian. "Flat" here
#' means no prior ADDED, not the fit's own penalty stripped - for the
#' bare likelihood of a MAP fit, write `frm_sample(fit, prior =
#' "flat")` instead, and read it as a different question.
#'
#' What this function never adds is [frm_sample()]'s DEFAULT priors.
#' Those defaults
#' apply on both of `frm_sample()`'s routes, and this function opts out
#' of them explicitly (`frm_sample(.diagnostic = TRUE)`) rather than
#' inheriting whatever the default is: a prior nothing in the fit ever
#' saw would change the very thing being measured.
#'
#' One density is added, on a model with `mo()` terms, and it adds
#' nothing to the comparison: each simplex is sampled in softmax
#' coordinates under `dirichlet(1)`, the flat density on the simplex, so
#' that the posterior of the weights is the likelihood the fit maximized
#' (times a MAP fit's own prior). With no density at all the sampled
#' coordinates would be the prior, which in softmax coordinates is
#' `1 / prod(w)` and in the fit's own coordinates is not integrable. The
#' simplex is compared on the weight scale: its rows are brms's `simo_`
#' weights, with `ml` and `wald_se` the weights and delta-method errors
#' `summary()` reports. A weight at 0 or 1 (a face or a vertex of the
#' simplex) has no Wald error (`NA`) and is reported in a message rather
#' than flagged.
#'
#' That is also why it samples CENTERED on an ordinary fit.
#' `frm_sample()`'s non-centered
#' parameterization (see Reparameterization there) is offered only for
#' blocks whose variance parameters carry a prior, and with the
#' defaults off none do: the
#' flat prior that makes the comparison meaningful is exactly the one
#' that leaves a flat tail at `sd = 0` for a non-centered chain to walk
#' into. So the reparameterization default costs this function nothing
#' and changes nothing about it. Give the variance parameters a prior
#' through `prior =`
#' and the run non-centers; but then it is measuring the Laplace
#' approximation of a different posterior, which is usually not the
#' question.
#'
#' @param fit A `frmtmb_fit`.
#' @param chains,iter Passed to [frm_sample()].
#' @param ... Passed to [frm_sample()].
#' @return A data frame (one row per outer parameter) with columns
#'   `ml`, `post_mean`, `wald_se`, `post_sd`, `z_shift`
#'   ((post_mean - ml)/post_sd) and `sd_ratio` (post_sd/wald_se).
#'
#' @details
#' The srr standard this function used to carry (RE1.4, "the
#' assumptions the model rests on are documented and testable") is
#' claimed by `frmtmb::frm()` now that the two live in different
#' packages: the assumption is the FIT's, and core has to state it
#' whether or not this package is installed. See the Laplace
#' approximation section of `?frmtmb::frm`.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#' # a binary GLMM with small clusters: the regime where the Laplace
#' # approximation and Wald intervals are least reliable
#' set.seed(4)
#' dd <- data.frame(x = rnorm(120), g = factor(rep(1:30, 4)))
#' dd$y <- rbinom(120, 1,
#'                plogis(0.3 + 0.5 * dd$x + rnorm(30, 0, 1)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + frmtmb::bernoulli(), data = dd)
#'
#' cl <- check_laplace(fit, chains = 1, iter = 500, refresh = 0)
#' cl
#' # |z_shift| well above 0 or sd_ratio far from 1 marks the parameters
#' # whose Wald interval to replace with a profile or bootstrap one
#' cl[abs(cl$z_shift) > 0.3 | cl$sd_ratio > 1.3, ]
#' }
#' }
#' @export
check_laplace <- function(fit, chains = 2, iter = 1000, ...) {
  # the whole point is NUTS against the ML mode and its Wald standard
  # errors; without a mode there is nothing to check the sampler against
  require_fitted(fit, "check_laplace()")
  # explicit, not inherited: frm_sample() defaults to the brms priors on
  # both of its routes now, and this function must keep sampling the
  # density the fit MAXIMIZED - otherwise it would compare NUTS on one
  # posterior against a mode and Wald errors from another. That density
  # is the flat likelihood for an ordinary fit and the fit's own
  # penalized one for a MAP fit, which is what the flag leaves standing
  ds <- frm_sample(fit, chains = chains, iter = iter, ...,
                   .diagnostic = TRUE)
  # the ML mode and its standard errors are on the link scale
  m <- draws_internal_matrix(ds)
  # the sd_ and cor_ columns are summaries of the theta beside them, not
  # parameters with an ML value of their own, and the internal matrix
  # has none of them
  keep <- setdiff(draws_outer_cols(ds),
                  draws_natural_cols(ds$fit)$derived$names)
  ml <- fit$opt$par
  se <- sqrt(diag(sdr_of(fit)$cov.fixed))
  # A mo() simplex is compared on the weight scale. Its coordinates are
  # not comparable: the fit's chart has 2^D sheets, and on 106 of 200
  # simplexes of ls ~ mo(income) * age (seeds 1 to 100) the fit sits on
  # another sheet than the one draws are carried back to (the review's
  # dev/optima-rev2-laplace.R). The ML weights and their delta-method
  # errors are summary()'s (summary_mo_frame()); the draws are the
  # simo_ columns. The coordinates leave both sides.
  mo_terms <- simplex_frame_terms(fit)
  mo_tab <- NULL
  if (length(mo_terms)) {
    in_mo <- names(ml) %in% vapply(mo_terms, `[[`, "", "zeta")
    ml <- ml[!in_mo]
    se <- se[!in_mo]
    mo_names <- unlist(lapply(mo_terms, `[[`, "names"))
    keep <- setdiff(keep, mo_names)
    mo_tab <- summary_mo_frame(fit, 0.95, mo_terms)
  }
  stopifnot(length(ml) == length(keep))
  post_mean <- colMeans(m[, keep, drop = FALSE])
  post_sd <- apply(m[, keep, drop = FALSE], 2, stats::sd)
  if (!is.null(mo_tab)) {
    W <- ds$draws[, mo_names, drop = FALSE]
    keep <- c(keep, mo_names)
    ml <- c(ml, mo_tab[, "Estimate"])
    se <- c(se, mo_tab[, "Est.Error"])
    post_mean <- c(post_mean, colMeans(W))
    post_sd <- c(post_sd, apply(W, 2, stats::sd))
    m <- cbind(m[, setdiff(colnames(m), mo_names), drop = FALSE], W)
  }
  # the check is only as good as the chain: a short chain that wandered
  # inflates post_sd and reads as "Laplace questionable" when the truth
  # is "chain unusable", so the effective sample size rides along and
  # an unhealthy chain is called out as such
  ess <- if (requireNamespace("posterior", quietly = TRUE)) {
    vapply(keep, function(nm) posterior::ess_bulk(m[, nm]), numeric(1))
  } else {
    rep(NA_real_, length(keep))
  }
  out <- data.frame(
    parameter = keep,
    ml = unname(ml),
    post_mean = unname(post_mean),
    wald_se = unname(se),
    post_sd = unname(post_sd),
    z_shift = unname((post_mean - ml) / post_sd),
    sd_ratio = unname(post_sd / se),
    ess_bulk = unname(ess),
    row.names = NULL
  )
  if (any(is.finite(out$ess_bulk) & out$ess_bulk < 100)) {
    frm_message("check_laplace(): the chain mixed too poorly to judge the ",
                "approximation (bulk ESS under 100 for ",
                paste(out$parameter[is.finite(out$ess_bulk) &
                                      out$ess_bulk < 100], collapse = ", "),
                "). Rerun with more iterations before reading z_shift or ",
                "sd_ratio")
  }
  flagged <- abs(out$z_shift) > 0.5 | out$sd_ratio > 1.5 |
    out$sd_ratio < 2 / 3
  # a weight at 0 or 1 has no Wald error (summary_mo_frame()), and no
  # Wald interval to judge: it is a boundary estimate, named as such
  # (its draws' name: simo_ for a mo() simplex, sbhaz for a baseline)
  face <- attr(mo_tab, "face_full") %||% character(0)
  flagged <- flagged %in% TRUE & !out$parameter %in% face
  if (any(flagged)) {
    frm_message("Laplace/Wald approximation questionable for: ",
                paste(out$parameter[flagged], collapse = ", "))
  }
  if (length(face)) {
    frm_message("check_laplace(): ", paste(face, collapse = ", "),
                " sit", if (length(face) == 1L) "s", " at 0 or 1, on the ",
                "boundary of the simplex, where the Laplace approximation ",
                "and the Wald error do not apply; read the posterior ",
                "instead")
  }
  out
}

#' Sample from a frmtmb fit with tmbstan (NUTS)
#'
#' Hands the fitted RTMB object to [tmbstan::tmbstan()]; from Stan's
#' point of view the model is the (Laplace-free) joint density, so all
#' parameters including random effects are sampled unless `laplace =
#' TRUE` is passed through.
#'
#' It hands over the fit's OWN objective and nothing else: no default
#' priors, no non-centering, no named draws. [frm_sample()] is the
#' route that applies brms's default priors and returns the draws
#' surface; this one is the escape hatch to tmbstan's own arguments.
#' A `mo()` simplex is therefore sampled flat in the fit's own
#' coordinates, with no density on it, and those coordinates fold at the
#' simplex's faces and have one point at infinity: on a weakly
#' identified simplex the draws collapse there (to the uniform simplex,
#' with coordinates near 1e18). The function says so in a message on a
#' fit with `mo()` terms; [frm_sample()] samples the simplex under
#' brms's `dirichlet(1)` instead. A MAP fit whose prior has a `simo` row
#' has that Dirichlet on its tape, with the Jacobian of its own
#' coordinates, so that simplex has a proper density here and the
#' message leaves it out; the coordinates still fold at the faces, which
#' keeps a chain on one of their equivalent sheets.
#' `control` means here what it means in [frm_sample()] and in brms,
#' the sampler's own list, because everything reaches
#' [tmbstan::tmbstan()] unchanged.
#'
#' Two things are NOT an escape hatch, and both are failures this
#' function used to pass on in silence.
#'
#' The first is the compatibility registry's refusals. This function
#' runs the same pre-flight [frm_sample()] does (see its "The pre-flight
#' against the compatibility registry" section), because a model the
#' sampler cannot run fails harder here, and the abort that produces the
#' failure can leave rstan unusable for the rest of the R session.
#'
#' The second is an empty result. rstan answers some failures by
#' declining to sample, printing its reason to the console and returning
#' a `stanfit` with nothing in it; this function returned that object.
#' It now stops instead, with the same message [frm_sample()] raises,
#' naming this function. The two refusals are deliberately different
#' sentences: the pre-flight one says the model cannot be sampled at
#' all, and this one says the run in front of you produced nothing.
#'
#' @param fit A `frmtmb_fit`.
#' @param ... Passed to [tmbstan::tmbstan()] (chains, iter, laplace,
#'   `control`, ...).
#' @return A `stanfit` object holding draws. **It is never an empty
#'   one:** rstan answers some failures by declining to sample, printing
#'   its reason and returning the shell of a `stanfit`, and until 0.3.2
#'   this function passed that shell back, so a failed run looked like a
#'   successful one until an accessor came back empty. It now stops with
#'   a message naming the run, which is a change to this function's
#'   return contract. A model that cannot be sampled at all is a
#'   different refusal and is raised earlier; see Details.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#' set.seed(9)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # The raw stanfit, for rstan and bayesplot code that wants one.
#' # Use frm_sample() instead when you want frmtmb parameter names.
#' # This run is deliberately short, so expect sampler warnings.
#' sf <- as_tmbstan(fit, chains = 1, iter = 400, refresh = 0)
#' class(sf)
#' # every parameter is sampled, random effects included, because Stan
#' # sees the joint density. Pass laplace = TRUE to integrate them out.
#' dim(as.matrix(sf))
#' }
#' }
#' @export
as_tmbstan <- function(fit, ...) {
  # The pre-flight belongs on BOTH doors, and this is the door where a
  # refused model does more damage: frm_sample() at least catches an
  # empty stanfit, while this function hands one back without a word,
  # and the abort that produced it leaves rstan's nested autodiff arena
  # unbalanced for the rest of the session ("empty_nested() must be
  # true before calling recover_memory()"), so every later chain in
  # that session fails too. Verified end to end on an frm_ode() fit.
  sample_preflight(fit, what = "as_tmbstan()")
  if (!requireNamespace("tmbstan", quietly = TRUE)) {
    frm_stop("as_tmbstan() needs the 'tmbstan' package", call. = FALSE)
  }
  check_tmbstan_build("as_tmbstan()")
  # a MAP fit's simo prior is on its tape, in the fit's own chart, and
  # gives that simplex a proper density there; only the others are flat
  mo <- simplex_frame_terms(fit)
  priored <- mo_fit_priored(fit)
  flat <- Filter(function(t) !t$zeta %in% priored, mo)
  if (length(flat)) {
    frm_message("as_tmbstan(): the simplex ",
                paste(sub("[[]1[]]$", "", vapply(flat, function(t) {
                  t$names[1L]
                }, "")), collapse = ", "),
                " has no density on this route and is sampled flat in ",
                "the fit's own coordinates, which a weakly identified ",
                "simplex collapses in. frm_sample() samples it under ",
                "brms's dirichlet(1)")
  }
  sf <- tmbstan::tmbstan(fit$obj, ...)
  check_stan_draws(sf, "as_tmbstan()")
  sf
}

#' Refuse an empty stanfit rather than hand one back.
#'
#' tmbstan can return a stanfit with no draws in it AND no error: rstan
#' declines to sample, prints its reason to the console, and gives back
#' the shell. `frm_sample()` has always refused that. `as_tmbstan()`
#' returned it, so the escape hatch answered a failed run with an object
#' whose every accessor is empty and whose emptiness the caller had to
#' notice for themselves. Both call this now, and the message names the
#' caller.
#'
#' The message says what the failure is NOT, because there are two ways
#' to get nothing out of this package and they want different fixes: a
#' model the sampler cannot run at all is refused by the compatibility
#' pre-flight before any chain starts, and this one is a failure of the
#' run in front of you.
#'
#' The cause is left generic on purpose. The known case is a tape that
#' calls an external solver, and the package supplying such a tape
#' documents its own failure; naming that package from here would tie
#' this message to something frmtmb.sample does not depend on, for a
#' pointer that package's own documentation gives better.
#'
#' On the reads: `@` cannot partial-match (R raises "no slot of name"),
#' and `[[` rather than `$` on `sim` is the house rule, not a fix.
#' Nothing was reachable through the old `$`: rstan's `sim` carries
#' twelve names and none but `samples` begins with "samples", and on a
#' failed run `sim` is empty, so `||` short-circuits and the second
#' test is never evaluated.
#'
#' @noRd
check_stan_draws <- function(sf, what) {
  if (length(sf@sim) && length(sf@sim[["samples"]])) {
    return(invisible(NULL))
  }
  frm_stop(what, ": the sampler returned no draws (rstan printed the cause ",
           "above). This is a failure of THIS RUN, not a model that cannot ",
           "be sampled at all: those are refused before any chain starts, ",
           "and frm_compat(\"frm_sample\") lists them. A known case here: ",
           "a tape that calls an external solver can fail inside tmbstan ",
           "even at the fitted optimum, where the same tape optimizes ",
           "without complaint; the package supplying such a tape documents ",
           "that. A solver failure mid-run also corrupts rstan's sampler ",
           "state, so retry in a fresh R session", call. = FALSE)
}

#' Refuse a tmbstan build that silently samples the wrong density.
#'
#' tmbstan generates its Stan model at INSTALL time: tools/autogen.R
#' rewrites a `std_normal_lpdf` placeholder in the stanc output into
#' the call that evaluates the TMB objective. stanc 2.39.0 (shipped by
#' StanHeaders 2.39.1, on CRAN 2026-09-02) emits TWO log_prob_impl
#' overloads where 2.32 emitted one, and autogen patches only the
#' first match. HMC reads value and gradient from the unpatched
#' reverse-mode overload, so every chain samples a standard normal
#' instead of the model, with the objective, priors and data all
#' silently absent. Proven by a one-line autogen patch flipping the
#' behavior; dev/prior-dropping-investigation.md has the full record.
#'
#' The check is static (one read of the installed model.hpp, cached
#' per session) and names the exact defect, so a user on an affected
#' build gets a refusal before any sampling rather than plausible
#' garbage after it.
#'
#' WHERE IT FAILS OPEN, which is deliberate and has to stay written
#' down. `nzchar(hpp) &&` makes an unreadable installation count as
#' healthy. Three shapes were measured through that door and all three
#' return FALSE, so the sampler runs (dev/tmbstan-findings.md):
#'
#' - no model.hpp at all;
#' - a model.hpp that is present but EMPTY;
#' - the defect written under a renamed placeholder, say
#'   `std_normal_lpdf<propto__, false>(y)`.
#'
#' Fail-open is right HERE and only here. Refusing every installation
#' whose model.hpp cannot be read would break correct machines, and
#' each reachable absence also breaks tmbstan itself: its `configure`
#' runs autogen at install time and `stopifnot(length(i) >= 1)` aborts
#' the build rather than shipping a package without the file, so the
#' user meets a loud runtime error and not plausible draws. A CALLER
#' THAT TURNS THIS ANSWER INTO A PRINTED CLAIM MUST NOT INHERIT IT:
#' the pin step in .github/workflows/check-frmtmb-sample.yaml keeps
#' absent and empty as outcomes of their own and fails the job on
#' them.
#'
#' A renamed placeholder is the shape no string check can cover. What
#' covers it is the entry point the defect actually replaces:
#' `rstan::grad_log_prob()` reads the unpatched reverse-mode overload,
#' so on an affected build it returns the standard normal's gradient
#' rather than the model's. tests/testthat/test-tmbstan-build-guard.R
#' asserts those two against each other and takes no position on any
#' string.
#'
#' @noRd
tmbstan_build_broken <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    hpp <- system.file("model.hpp", package = "tmbstan")
    cached <<- nzchar(hpp) &&
      any(grepl("std_normal_lpdf<propto__>(y)",
                readLines(hpp, warn = FALSE), fixed = TRUE))
    cached
  }
})

check_tmbstan_build <- function(caller) {
  if (tmbstan_build_broken()) {
    frm_stop(caller, ": this tmbstan installation was built against ",
             "StanHeaders >= 2.39, whose code generator leaves one of the ",
             "two generated log-density overloads unpatched, so EVERY ",
             "chain silently samples a standard normal instead of the ",
             "model (tmbstan tools/autogen.R replaces only the first ",
             "match). Reinstall tmbstan (>= 1.2.1), which patches every ",
             "overload, and distrust any draws already produced by this ",
             "installation",
             call. = FALSE)
  }
  invisible(NULL)
}
