# Why this file exists.
#
# frmtmb ports the brms vocabulary, so it exports 30 S3 generics, and
# 28 of those names are already owned by another package. An exported
# generic is not a neutral alias: UseMethod() consults the method table
# of the namespace where the generic IT REACHED was defined. A frmtmb
# generic sitting above brms on the search path therefore sends a
# `brmsfit` into frmtmb's table, which has no entry for that class, and
# the call stops at "no applicable method" while brms's own method sits
# registered and intact in brms's table. Measured on the base commit:
# after `library(brms); library(frmtmb)`, ALL 27 of the generics the
# two share lose dispatch, and two of the 27 lose it SILENTLY, falling
# into frmtmb's own `.default` rather than erroring.
#
# Delayed registration alone does NOT fix this, and that is the reason
# the defect survived so long: the base commit already carried 38
# `S3method(pkg::generic, class)` directives over 36 distinct generics,
# and the defect reproduces there. Putting frmtmb's method in the
# OWNER's table helps the OWNER's generic. It does nothing about
# frmtmb's rival one.
#
# So frmtmb must stop owning a name someone else owns. Two routes:
#
#  * Static. `fixef`, `ranef` and `VarCorr` come from nlme, which is a
#    Recommended package that ships with R, and `refit` from generics,
#    which declares only base packages. Both are importFrom()-ed and
#    re-exported, which is exactly what lme4 and brms do with nlme.
#    nlme is then the ONE generic that frmtmb, lme4, glmmTMB and brms
#    all dispatch through, because all three of the others import it
#    from nlme too.
#
#    `refit` is the exception in that group and needs both routes.
#    lme4 DEFINES its own `refit` rather than importing generics', so
#    `refit.merMod` lives in lme4's table and the generics import does
#    not reach it: measured, `refit` was the one name still lost to an
#    lme4 user after the import landed. So lme4 is in the run-time
#    table for `refit` as well, with generics' generic as the fallback.
#
#  * An ACTIVE BINDING, for the owners frmtmb cannot afford in Imports:
#    bayesplot pulls ggplot2 and brms pulls rstan, and the whole point
#    of frmtmb is to need neither.
#
# Why an active binding rather than swapping the binding at load.
#
# The first version of this file swapped frmtmb's binding for the
# owner's inside `.onLoad`, and set a load hook per owner for an owner
# that turned up later. It worked, and review broke it three ways that
# an active binding cannot be broken:
#
#  * The hook path had to write into a SEALED namespace and into the
#    attached `package:frmtmb`, so it called `unlockBinding()`, which
#    R CMD check reports as a possibly unsafe call. `makeActiveBinding()`
#    in `.onLoad` runs while the namespace is still unsealed, so there
#    is nothing to unlock and no NOTE.
#  * A binding swapped once is a hard reference to one namespace
#    INSTANCE. `library(loo); library(frmtmb)`, then unload loo, then
#    load anything that pulls a fresh loo, and frmtmb was dispatching
#    into a dead namespace's method table with nothing said. An active
#    binding re-resolves on every access, so a stale instance cannot
#    survive one call.
#  * Active-ness PROPAGATES through `importIntoEnv()`. The attached
#    `package:frmtmb` binding is active, and so is the binding in a
#    package that re-exports the name from frmtmb, which is how
#    `frmtmb.sample` gets the repair too. A swapped binding reached
#    neither, because both are copies taken at attach time.
#
# The cost is one closure call per ACCESS of the name instead of one
# per session, and it is a real cost rather than a rounding error:
# see dev/generics-adoptcost.R for the number on this box, measured
# against an ordinary binding lookup and against one `fixef()` call.
# These names are reached once per user action, not in a loop, so
# the trade is the right one; it is stated rather than waved away
# because the first spelling of this file called it small and the
# measurement did not agree.

#' Generics frmtmb re-exports from the package that owns them.
#'
#' nlme is Recommended and ships with R; generics declares only base
#' packages, `library(generics)` costing 0.000 s at the minimum of 12
#' replicates against 0.860 s for lme4 (`dev/generics-loadcost.R`).
#' So these four are free. Three of them need nothing further;
#' `refit` also joins the run-time table, because lme4 defines a
#' rival `refit` that this import cannot reach.
#'
#' @importFrom nlme fixef ranef VarCorr
#' @importFrom generics refit
#' @rawNamespace export(VarCorr,fixef,ranef,refit)
#' @name frmtmb-shared-generics
#' @keywords internal
NULL

#' The generic names frmtmb exports that another package owns, with the
#' owner to defer to.
#'
#' Ownership was audited with `parseNamespaceFile()` over every
#' installed package, never with grep, because nlme declares its
#' exports in a multi-name block that a one-name-per-line grep cannot
#' see. `dev/generics-audit.R` and `dev/generics-audit2.R` hold the
#' audit; the second also reports where each brms method is registered,
#' which is what decides the target: brms imports `loo` from loo and
#' `as_draws_df` from posterior, so its methods for those live in the
#' OWNER's table and deferring to the owner reaches them.
#'
#' `ngrps` and `refit` are the two names with more than one owner, and
#' the owners' generics are different closures whose methods sit only
#' in their own tables. They collide with each other with or without
#' frmtmb, so frmtmb picks whichever one the user would have reached
#' without it: the search path first, then this order.
#'
#' @noRd
frm_generic_owners <- list(
  as_draws = "posterior",
  as_draws_array = "posterior",
  as_draws_df = "posterior",
  as_draws_list = "posterior",
  as_draws_matrix = "posterior",
  as_draws_rvars = "posterior",
  nchains = "posterior",
  ndraws = "posterior",
  niterations = "posterior",
  nvariables = "posterior",
  variables = "posterior",
  loo = "loo",
  loo_compare = "loo",
  waic = "loo",
  bayes_R2 = "rstantools",
  prior_summary = "rstantools",
  pp_check = "bayesplot",
  conditional_effects = "brms",
  expose_functions = "brms",
  hypothesis = "brms",
  posterior_summary = "brms",
  LOO = "brms",
  WAIC = "brms",
  ngrps = c("brms", "lme4"),
  refit = "lme4"
)

#' Is this the owner's GENERIC, or only a function with the right name?
#'
#' A package that exported a plain function under one of these names
#' would otherwise take the binding, and frmtmb's method for its own
#' class would become unreachable with nothing said. All 26 owner and
#' name pairs are generics today; this is the absent case built rather
#' than assumed.
#'
#' @noRd
frm_is_generic <- function(f) {
  if (!is.function(f)) return(FALSE)
  b <- body(f)
  if (is.null(b)) return(FALSE)
  any(grepl("UseMethod", deparse(b), fixed = TRUE))
}

#' Which owner of `gen` would the user have reached without frmtmb?
#'
#' Attachment order decides that, so the search path is consulted
#' first, and the owner nearest the front of it wins. The `owners`
#' order is the tie-break for owners that are loaded but not
#' attached.
#'
#' Only `ngrps` and `refit` reach this, and it runs on every ACCESS
#' of those two names, so it is two vector operations rather than a
#' loop with a `sub()` in it: the loop spelling measured 49.13 us
#' per access against 1.70 us for a single-owner name.
#'
#' @noRd
frm_adopt_target <- function(gen, owners) {
  pos <- match(paste0("package:", owners), search())
  if (any(!is.na(pos))) {
    p <- owners[which.min(ifelse(is.na(pos), Inf, pos))]
    e <- as.environment(paste0("package:", p))
    if (exists(gen, envir = e, inherits = FALSE)) return(p)
  }
  for (p in owners) {
    if (!isNamespaceLoaded(p)) next
    g <- tryCatch(getExportedValue(p, gen), error = function(e) NULL)
    if (is.function(g)) return(p)
  }
  NULL
}

#' Install one active binding.
#'
#' A function of its own for two reasons. The loop variable is
#' captured by value rather than by reference, and the memo lives in
#' THIS closure rather than in a shared environment, so resolving a
#' name costs no key to build and no environment to search. The memo
#' holds the owner's NAMESPACE ENVIRONMENT beside the function, and
#' validating it is an environment identity test: a namespace that
#' was unloaded and reloaded is a different environment, so a stale
#' instance cannot be served. That is the failure the round-1
#' binding swap had, and it is closed here by construction rather
#' than by a hook.
#'
#' @noRd
frm_bind_generic <- function(ns, gen, owners, fallback) {
  force(gen)
  force(fallback)
  single <- length(owners) == 1L
  one <- if (single) owners[[1L]] else NULL
  memo_ns <- NULL
  memo_fn <- NULL
  makeActiveBinding(gen, function() {
    own <- if (single) {
      if (isNamespaceLoaded(one)) one else NULL
    } else {
      frm_adopt_target(gen, owners)
    }
    if (is.null(own)) return(fallback)
    # no tryCatch on this line: it runs on every access, a
    # handler stack costs about 7 us here, and asNamespace()
    # cannot fail behind the isNamespaceLoaded() above
    ons <- asNamespace(own)
    if (!identical(ons, memo_ns)) {
      g <- tryCatch(getExportedValue(own, gen),
                    error = function(e) NULL)
      # a name an owner exports need not be a GENERIC; if it is not,
      # frmtmb keeps its own or its method for its own class goes
      # unreachable with nothing said
      if (!frm_is_generic(g)) g <- NULL
      memo_ns <<- ons
      memo_fn <<- g
    }
    if (is.null(memo_fn)) fallback else memo_fn
  }, ns)
  invisible(TRUE)
}

#' Replace every borrowed generic in this namespace with an active
#' binding that resolves to the owner's generic when the owner is
#' there.
#'
#' Called from `.onLoad`, which is the only place it can be called.
#' `loadNamespace()` runs `.onLoad` at line 367, after
#' `registerS3methods()` at line 345 and before the namespace is
#' sealed. So every `S3method()` directive has already resolved
#' against frmtmb's own generic and landed in frmtmb's own table,
#' which is exactly where the fallback needs them, and there is still
#' no lock to defeat.
#'
#' The table is an argument because an extension has the same defect
#' for the generics IT defines, and must be repaired from its own
#' `.onLoad`, the one moment ITS namespace is unsealed. One
#' implementation with two callers rather than a copy per package, so
#' that a repair to the binding reaches both. Exported on
#' `?frmtmb-sampling-api`.
#'
#' @noRd
frm_install_generics <- function(pkgname = "frmtmb",
                                 owners = frm_generic_owners) {
  # Refused rather than skipped: a malformed table would install
  # nothing and say nothing, which is the defect back again.
  ok <- is.list(owners) && length(owners) > 0L &&
    !is.null(names(owners)) && all(nzchar(names(owners))) &&
    all(vapply(owners, function(o) {
      is.character(o) && length(o) >= 1L && !anyNA(o) && all(nzchar(o))
    }, NA))
  if (!ok) {
    frm_stop("frm_install_generics(owners =) must be a named list of ",
             "non-empty character vectors: one entry per generic name, ",
             "naming the packages that own it in the order to prefer",
             call. = FALSE)
  }
  ns <- asNamespace(pkgname)
  done <- character()
  for (gen in names(owners)) {
    # Idempotent. A second call must not capture the LIVE generic as
    # the fallback, which would pin whatever happened to be loaded.
    if (exists(gen, envir = ns, inherits = FALSE) &&
        bindingIsActive(gen, ns)) {
      done <- c(done, gen)
      next
    }
    fb <- tryCatch(get(gen, envir = ns), error = function(e) NULL)
    if (!is.function(fb)) next
    # `refit` arrives through importFrom(), so it has no binding of its
    # own here and there is nothing to remove; the active binding
    # shadows the import instead.
    if (exists(gen, envir = ns, inherits = FALSE)) {
      rm(list = gen, envir = ns)
    }
    frm_bind_generic(ns, gen, owners[[gen]], fb)
    done <- c(done, gen)
  }
  invisible(done)
}
