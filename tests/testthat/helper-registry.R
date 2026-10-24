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
# The helper is frmtmb's exported frmtmb_local_registries() (R/
# registries.R), so that an extension's suite saves and restores the
# same registries this suite does. It covers every registry at once, so
# a test need not know which ones a registration touches, and after the
# restore it runs the .onLoad() of every frmtmb.* namespace that loaded
# inside the scope, in the order they loaded: the restore removed their
# entries, and a contributed rule's precedence is its position.
local_frmtmb_registries <- function(env = parent.frame()) {
  frmtmb::frmtmb_local_registries(env)
}

# The registries the helper saves, and the environments it leaves alone
# because they are caches. test-registry-helper.R asserts that the two
# lists are every top-level environment of the namespace, so a new or
# renamed registry fails a test until R/registries.R covers it.
frmtmb_test_registries <- frmtmb:::frmtmb_registry_envs$registries
frmtmb_test_caches <- frmtmb:::frmtmb_registry_envs$caches

# One registry environment, read so that a renamed one stops the test:
# `ns$<absent>` is NULL, and assigning into NULL modifies a local copy
# without a word.
frmtmb_test_registry <- function(nm) {
  e <- get(nm, envir = asNamespace("frmtmb"), inherits = FALSE)
  stopifnot(is.environment(e))
  e
}
