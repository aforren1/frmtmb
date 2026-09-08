# The hazard containers, and the scanner that keeps `$` off them.
#
# `$` on a list falls back to PARTIAL matching when no name matches
# exactly, so a slot that is absent does not read `NULL`, it reads a
# neighbor whose name it is a prefix of. The containers below are keyed
# by an open or prefix-colliding vocabulary, and the package has been
# bitten three times: `pars$b` read `beta`, `ctx$mix` read `mix_g`, and
# at v0.49 `par_template$theta` read `thetaac` on a model with residual
# autocorrelation and no random effects, so get_prior() offered rows for
# parameters the model did not have. `[[` has no fallback and returns
# NULL, which is what every one of those call sites already assumed it
# was getting.
#
# The list lives HERE, in R/, rather than in the test that reads it,
# because the same list has to reach every extension. It is the whole
# reason frm_hazard_reads() is exported: seven extension test suites
# calling one scanner cannot drift, and seven copies of the list can.
# dev/bracket-sweep.md holds the inventory it was drawn from.

#' The container names `$` is refused on. Not a control: adding a name
#' reserves it across the monorepo, and the sweep that placed these
#' collected the names each container is ever accessed with rather than
#' assuming them.
#'
#' WHAT IS NOT LISTED. Environments are absent on purpose: `$` on an
#' environment is already an exact lookup with no partial fallback, so
#' `frmtmb_aterm_registry`, `frmtmb_frame_checks` and `hyp_shadow_state`
#' are safe as written. `covstruct_registry` is a plain list and is
#' listed. Structs with a closed slot vocabulary (`resp`, `rspec`,
#' `spec`) and the user-facing S3 objects (`fit`, `x`, `object`) are not
#' listed either: their names are fixed, exact and documented.
#'
#' @noRd
hazard_containers <- c(
  # parameter and estimate lists: theta/thetaac/thetar, b/beta/betad
  "pars", "tpl", "par_template", "est", "estimates",
  # aterm containers: cens/cens_y2, mi/mi_sd, se/se_sigma
  "aterms", "av",
  # random-effects blocks: aux_D/aux_D2, aux_Q/aux_Qk
  "blk", "bk",
  # structure blocks, and the linear-predictor, autocor and eta-design
  # blocks: par/param_colnames, smooth/smooths, p/patterns, n/nonest
  "block", "blocks", "lp", "ac", "ed",
  # the frame itself: y/y_levels, linpred/linpreds
  "frame",
  # family objects: family/family_finalize, mix/mix_groups, sim/sim_ctx
  "fam", "family",
  # simulator contexts, and the distributional-parameter containers,
  # both keyed by names the user supplies
  "ctx", "dpars", "dpv", "dp",
  # registries that are lists rather than environments: us/us_t
  "covstruct_registry"
)

#' Find `$` reads of a partial-matching container
#'
#' A testing aid for frmtmb and its extensions, not part of a fitting
#' workflow. It reports every place in a package where `$` is used on
#' one of the containers whose slot names collide under partial
#' matching, so an extension can assert in its own test suite that it
#' has none.
#'
#' @details
#' `$` on a list partial-matches: when no name matches exactly, R
#' returns the slot whose name the requested one is a prefix of, rather
#' than `NULL`. frmtmb keeps several containers whose vocabulary is open
#' (a distributional parameter is named by the user) or prefix-colliding
#' by construction (`theta` against `thetaac`, `b` against `beta`, `se`
#' against `se_sigma`), and a `$` read of an absent slot on one of those
#' silently returns a neighbor. `[[` returns `NULL` instead, which is
#' what such call sites assume they are getting.
#'
#' The rule is keyed on the container's NAME, not on `$` in general:
#' `fit$obj`, `x$fit` and `object$draws` are untouched, and a chain is
#' judged one link at a time by its own left-hand name, which is why
#' `fit$frame[["par_template"]]` is right and `fit$frame$par_template`
#' is not. The names are therefore RESERVED: binding one of them to
#' something that is not the container it names is itself a hit, and the
#' fix is to rename the local rather than to exempt it.
#'
#' The scan reads the abstract syntax tree of every function in the
#' namespace, so it sees `$` in code and not in a string or a comment,
#' and it works against an installed package, where the sources are no
#' longer on disk. What it cannot see is top-level package code that is
#' not inside a function, because only that code's result is installed.
#'
#' @param package A package name, or an environment holding functions
#'   (a namespace, which is what a package name resolves to).
#' @return A character vector, one entry per read, spelled
#'   `"function: container$slot"`, sorted by function name. A package
#'   with no hazard read gives `character(0)`.
#' @examples
#' # in a package's own test suite:
#' #   expect_identical(frmtmb::frm_hazard_reads("frmtmb.spline"),
#' #                    character(0))
#' frm_hazard_reads("frmtmb")
#' @export
frm_hazard_reads <- function(package) {
  ns <- if (is.environment(package)) package else {
    if (!is.character(package) || length(package) != 1L || is.na(package)) {
      stop("`package` must be one package name, or an environment of ",
           "functions to scan", call. = FALSE)
    }
    getNamespace(package)
  }
  out <- character(0)
  # An empty argument (`x[, 1]`) errors when TOUCHED rather than when
  # extracted, so every inspection of a call's parts is guarded.
  #
  # A `function` written INSIDE a body carries its formals as a
  # pairlist, which is not a call, so a walk that descends only into
  # calls steps straight over an inner function's default arguments.
  # A default argument is code in whichever frame it runs, and the
  # source scanner reads it, so both scanners have to see it.
  part <- function(e, i) tryCatch({
    a <- e[[i]]
    if (is.pairlist(a) && length(a)) as.call(c(quote(list), as.list(a)))
    else if (is.call(a)) a
  }, error = function(err) NULL)
  # What the source-tree scanner calls "the token immediately to the
  # left of the `$`", as a name. A chain is judged one link at a time:
  # in `fit$frame$data` the left token of the second `$` is `frame`, so
  # that read is a hit even though the outer expression's own left side
  # is a call. `fit[["frame"]]$data` is not, and neither is `f(x)$data`,
  # because their left token is a bracket.
  lhs <- function(x) {
    if (is.name(x)) return(as.character(x))
    if (is.call(x) && identical(x[[1L]], quote(`$`)) && length(x) >= 3L &&
        is.name(x[[3L]])) {
      return(as.character(x[[3L]]))
    }
    NA_character_
  }
  walk <- function(e, where) {
    if (identical(e[[1L]], quote(`$`)) && length(e) >= 3L &&
        isTRUE(lhs(e[[2L]]) %in% hazard_containers)) {
      out[[length(out) + 1L]] <<- paste0(
        where, ": ", lhs(e[[2L]]), "$",
        paste(deparse(e[[3L]]), collapse = ""))
    }
    for (i in seq_along(e)) {
      a <- part(e, i)
      if (!is.null(a)) walk(a, where)
    }
    invisible(NULL)
  }
  for (nm in sort(ls(ns, all.names = TRUE))) {
    f <- tryCatch(get(nm, envir = ns, inherits = FALSE),
                  error = function(err) NULL)
    if (!is.function(f)) next
    b <- body(f)
    if (is.call(b)) walk(b, nm)
    # a default argument is code too, and runs in the callee's frame
    for (d in as.list(formals(f))) {
      if (tryCatch(is.call(d), error = function(err) FALSE)) walk(d, nm)
    }
  }
  out
}
