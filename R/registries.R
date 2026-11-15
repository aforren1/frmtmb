# The registries a package fills from its .onLoad(), saved and restored
# around a test, exported so that an extension's suite uses the same
# helper as frmtmb's own (tests/testthat/helper-registry.R wraps it).

#' The environments frmtmb keeps registrations in, and the ones that are
#' caches and are left alone: compat_accepts_cache is validated against
#' the contributed features, and compat_hand_cache and nl_shadow_cache
#' hold nothing a registration writes. test-registry-helper.R asserts
#' that the two lists are every top-level environment of the namespace,
#' so a new or renamed registry fails a test until it is listed here.
#'
#' @noRd
frmtmb_registry_envs <- list(
  registries = c("frmtmb_compat_contrib", "frmtmb_compat_cache",
                 "frmtmb_aterm_registry", "frmtmb_frame_checks",
                 "frmtmb_prior_defaults", "frmtmb_notice_state",
                 "frmtmb_covstruct_registry",
                 "frmtmb_autocor_registry"),
  caches = c("compat_accepts_cache", "compat_hand_cache", "nl_shadow_cache"))

#' `pkgs` in an order that puts a package after the frmtmb packages of
#' `pkgs` it imports (frmtmb.learn after frmtmb.eam).
#'
#' @noRd
registry_load_order <- function(pkgs) {
  out <- character()
  while (length(pkgs)) {
    ready <- vapply(pkgs, function(p) {
      !any(names(getNamespaceImports(p)) %in% setdiff(pkgs, p))
    }, TRUE)
    if (!any(ready)) ready[[1L]] <- TRUE
    out <- c(out, pkgs[ready])
    pkgs <- pkgs[!ready]
  }
  out
}

#' The key a registered function is known by, or a fresh environment
#' when its sameness cannot be proved. A `pkgload::load_all()` of an
#' extension runs its `.onLoad()` again with new closures, so identity
#' cannot be the key; two closures one factory made for different values
#' do different things, so code alone cannot be it either. A merge drops
#' a check with no word and a duplicate only runs one twice, so the key
#' merges only what is provably the same: the same code, and the same
#' values in every environment from the function's own up to the first
#' namespace, the global environment, base or empty. An environment or
#' an external pointer among those values is written by `serialize()` by
#' its contents or not at all, so equal bytes would not prove the same
#' object; meeting one makes the key unique. So does meeting any other
#' named environment on the walk (`attr(e, "name")`, an `attach()`ed
#' one): its name says nothing of its values, and above an attached one
#' lies the search path.
#'
#' @noRd
registration_key <- function(fn) {
  attr(fn, "frmtmb_registration_key") <- NULL
  env <- environment(fn) %||% baseenv()
  owner <- environmentName(topenv(env))
  seen <- new.env(parent = emptyenv())
  seen$ref <- FALSE
  hook <- function(x) {
    seen$ref <- TRUE
    NULL
  }
  frames <- list()
  while (!registration_named_env(env)) {
    if (nzchar(environmentName(env))) return(new.env(parent = emptyenv()))
    frames[[length(frames) + 1L]] <- tryCatch(
      serialize(as.list(env, all.names = TRUE, sorted = TRUE), NULL,
                refhook = hook),
      error = function(e) {
        seen$ref <- TRUE
        NULL
      })
    if (seen$ref) return(new.env(parent = emptyenv()))
    env <- parent.env(env)
  }
  # deparse() leaves the srcref out by default, so a load_all() closure
  # and an installed one of the same code share a key
  list(owner = owner, code = paste(deparse(fn), collapse = "\n"),
       frames = frames)
}

#' Where the walk of `registration_key()` stops: an environment a reload
#' does not copy, so that the same closure reloaded meets the same one.
#'
#' @noRd
registration_named_env <- function(env) {
  isNamespace(env) || identical(env, globalenv()) ||
    identical(env, baseenv()) || identical(env, emptyenv())
}

#' `fns` with `fn` registered once: a function whose key matches an
#' entry's, or the very same function, replaces that entry in place, so a
#' reloaded extension's check runs once per fit and keeps its position,
#' which is its order. Anything else is appended. Each entry carries its
#' key, computed once, so a registration does not serialize every other
#' entry's captured values again.
#'
#' @noRd
register_once <- function(fns, fn) {
  key <- registration_key(fn)
  plain <- fn
  attr(plain, "frmtmb_registration_key") <- NULL
  attr(fn, "frmtmb_registration_key") <- key
  for (i in seq_along(fns)) {
    old <- fns[[i]]
    k <- attr(old, "frmtmb_registration_key", exact = TRUE)
    attr(old, "frmtmb_registration_key") <- NULL
    if (identical(old, plain) || (!is.null(k) && identical(k, key))) {
      fns[[i]] <- fn
      return(fns)
    }
  }
  c(fns, list(fn))
}

#' Save frmtmb's registries, and restore them when a test ends
#'
#' A test that registers a covariance structure, a residual correlation
#' term, an addition term, a family's compatibility rows, a frame check
#' or a prior default calls `frmtmb_local_registries()` first. When
#' `env` exits (the end of the `test_that()` block), every registry is
#' put back as it was, so the registration cannot reach the tests after
#' it: `R CMD check` and covr run every file of a suite in one session.
#'
#' The restore puts back what it found. It never deletes an entry by
#' name, because the name may belong to a package that was loaded
#' before the test. A frmtmb package that loads inside the scope keeps
#' its registrations: after the restore, the `.onLoad()` of every
#' `frmtmb.*` namespace that appeared in the scope runs again, in the
#' order they loaded. A package load hook ([setHook()]) that other code
#' sets inside the scope stays set.
#'
#' @param env The environment whose exit restores the registries; the
#'   default is the caller's, which in a `test_that()` block is the test.
#' @return The saved state, invisibly. Called for the restore it
#'   schedules.
#' @seealso [frmtmb_register_covstruct()], [frmtmb_register_autocor()],
#'   [frmtmb_register_aterm()], [frmtmb_register_compat()] and
#'   [frmtmb_register_frame_check()], the registrations it undoes.
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' f <- function() {
#'   frmtmb_local_registries()
#'   # a registration here is undone when f() returns
#'   invisible(NULL)
#' }
#' f()
#' @export
frmtmb_local_registries <- function(env = parent.frame()) {
  if (!requireNamespace("withr", quietly = TRUE)) {
    frm_stop("frmtmb_local_registries() needs the withr package, which ",
             "testthat installs", call. = FALSE)
  }
  contrib <- frmtmb_compat_contrib
  cache <- frmtmb_compat_cache
  atreg <- frmtmb_aterm_registry
  checks <- frmtmb_frame_checks
  priors <- frmtmb_prior_defaults
  notices <- frmtmb_notice_state
  covs <- frmtmb_covstruct_registry
  acs <- frmtmb_autocor_registry
  saved <- list(
    features = contrib$features, rules = contrib$rules,
    cache_features = cache$features, cache_tbl = cache$tbl,
    atreg = atreg$reg, frame_checks = checks$fns,
    priors = priors$providers, covstructs = covs$reg,
    autocors = acs$reg,
    notices = as.list(notices, all.names = TRUE),
    loaded = loadedNamespaces())
  # record the order the frmtmb packages not yet loaded load in; the
  # hooks they had are put back when the scope ends
  seen <- new.env(parent = emptyenv())
  seen$order <- character()
  cand <- setdiff(unique(list.files(.libPaths(), "^frmtmb[.]")),
                  saved$loaded)
  mine <- lapply(cand, function(p) {
    force(p)
    function(...) seen$order <- c(seen$order, p)
  })
  # the restore removes only this scope's own hooks, by identity: a hook
  # another party set inside the scope is theirs and has to survive it.
  # It is scheduled before any hook is set, so an error while setting
  # them cannot leave one behind.
  withr::defer({
    for (i in seq_along(cand)) {
      ev <- packageEvent(cand[[i]], "onLoad")
      now <- getHook(ev)
      keep <- Filter(function(h) !identical(h, mine[[i]]), now)
      setHook(ev, keep, "replace")
    }
    contrib$features <- saved$features
    contrib$rules <- saved$rules
    cache$features <- saved$cache_features
    cache$tbl <- saved$cache_tbl
    atreg$reg <- saved$atreg
    checks$fns <- saved$frame_checks
    priors$providers <- saved$priors
    covs$reg <- saved$covstructs
    acs$reg <- saved$autocors
    rm(list = ls(notices, all.names = TRUE), envir = notices)
    list2env(saved$notices, envir = notices)
    new <- grep("^frmtmb[.]", setdiff(loadedNamespaces(), saved$loaded),
                value = TRUE)
    # in load order; a namespace no hook saw (a package outside
    # .libPaths()) goes after them, after its frmtmb imports
    run <- c(intersect(seen$order, new),
             registry_load_order(setdiff(new, seen$order)))
    for (p in run) {
      hook <- get0(".onLoad", envir = asNamespace(p), inherits = FALSE)
      if (is.function(hook)) hook(dirname(getNamespaceInfo(p, "path")), p)
    }
  }, envir = env)
  for (i in seq_along(cand)) {
    setHook(packageEvent(cand[[i]], "onLoad"), mine[[i]])
  }
  invisible(saved)
}
