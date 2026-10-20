# Save frmtmb's session registries and restore them when `env` exits.
#
# covr and R CMD check run every test file in ONE session, so a
# registration a test makes, or an entry it deletes, reaches every file
# after it. 0.71.0's coverage job failed that way: test-parse.R removed
# the dec() entry frmtmb.eam had registered when an earlier file loaded
# it. A test that registers saves first and puts back what it found,
# never deleting by name, because the name may belong to an extension
# that was loaded before the test ran.
#
# Every registry is covered at once, so a test need not know which ones
# a registration touches: frmtmb_register_aterm() also writes a compat
# feature, and the vocabulary cache follows the features. The notice
# state is here because a test that arms a notice clears it.
#
# An extension that LOADS inside the scope registers from its .onLoad(),
# and the restore would remove those entries while the namespace stays
# loaded: frmtmb.eam loaded and dec() unregistered, which is 0.71.0's
# failure the other way round. So after the restore the helper runs the
# registration again: it calls the .onLoad() of every frmtmb.* namespace
# that appeared in the scope, in the order they loaded. The order is
# recorded by an onLoad hook the helper sets for each frmtmb package not
# yet loaded, because a contributed rule's precedence is its position
# and loadedNamespaces() does not keep load order (a rerun in another
# order gave a different table). The restore removed their entries, so
# each lands once. Every extension's .onLoad() only registers, apart
# from frmtmb.sample's frm_install_generics(), which is idempotent.

# The registries the helper saves, and the environments it leaves alone
# because they are caches: compat_accepts_cache is validated against the
# contributed features, and compat_hand_cache and nl_shadow_cache hold
# nothing a test can register. test-registry-helper.R asserts that the
# two lists are every top-level environment of the namespace, so a new
# or renamed registry fails a test until it is covered here.
frmtmb_test_registries <- c(
  "frmtmb_compat_contrib", "frmtmb_compat_cache", "frmtmb_aterm_registry",
  "frmtmb_frame_checks", "frmtmb_prior_defaults", "frmtmb_notice_state")
frmtmb_test_caches <- c("compat_accepts_cache", "compat_hand_cache",
                        "nl_shadow_cache")

# One registry environment, read so that a renamed one stops the test:
# `ns$<absent>` is NULL, and assigning into NULL modifies a local copy
# without a word.
frmtmb_test_registry <- function(nm) {
  e <- get(nm, envir = asNamespace("frmtmb"), inherits = FALSE)
  stopifnot(is.environment(e))
  e
}

# `pkgs` in an order that puts a package after the frmtmb packages of
# `pkgs` it imports (frmtmb.learn after frmtmb.eam).
frmtmb_test_load_order <- function(pkgs) {
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

local_frmtmb_registries <- function(env = parent.frame()) {
  contrib <- frmtmb_test_registry("frmtmb_compat_contrib")
  cache <- frmtmb_test_registry("frmtmb_compat_cache")
  aterms <- frmtmb_test_registry("frmtmb_aterm_registry")
  checks <- frmtmb_test_registry("frmtmb_frame_checks")
  priors <- frmtmb_test_registry("frmtmb_prior_defaults")
  notices <- frmtmb_test_registry("frmtmb_notice_state")
  saved <- list(
    features = contrib$features, rules = contrib$rules,
    cache_features = cache$features, cache_tbl = cache$tbl,
    aterms = aterms$reg, frame_checks = checks$fns,
    priors = priors$providers,
    notices = as.list(notices, all.names = TRUE),
    loaded = loadedNamespaces())
  # record the order the frmtmb packages not yet loaded load in; the
  # hooks they had are put back when the scope ends
  seen <- new.env(parent = emptyenv())
  seen$order <- character()
  cand <- setdiff(unique(list.files(.libPaths(), "^frmtmb[.]")),
                  saved$loaded)
  old_hooks <- lapply(cand, function(p) getHook(packageEvent(p, "onLoad")))
  for (p in cand) {
    setHook(packageEvent(p, "onLoad"), local({
      q <- p
      function(...) seen$order <- c(seen$order, q)
    }))
  }
  withr::defer({
    for (i in seq_along(cand)) {
      setHook(packageEvent(cand[[i]], "onLoad"), old_hooks[[i]], "replace")
    }
    contrib$features <- saved$features
    contrib$rules <- saved$rules
    cache$features <- saved$cache_features
    cache$tbl <- saved$cache_tbl
    aterms$reg <- saved$aterms
    checks$fns <- saved$frame_checks
    priors$providers <- saved$priors
    rm(list = ls(notices, all.names = TRUE), envir = notices)
    list2env(saved$notices, envir = notices)
    new <- grep("^frmtmb[.]", setdiff(loadedNamespaces(), saved$loaded),
                value = TRUE)
    # in load order; a namespace no hook saw (a package outside
    # .libPaths()) goes after them, after its frmtmb imports
    run <- c(intersect(seen$order, new),
             frmtmb_test_load_order(setdiff(new, seen$order)))
    for (p in run) {
      hook <- get0(".onLoad", envir = asNamespace(p), inherits = FALSE)
      if (is.function(hook)) hook(dirname(getNamespaceInfo(p, "path")), p)
    }
  }, envir = env)
  invisible(saved)
}
