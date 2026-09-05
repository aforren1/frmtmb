# What frmtmb_register_compat() refuses.
#
# The registry answers for features that live in other packages, and it
# can only do that if the contributed rules reach the pairs they are
# written for. A rule side that names nothing matches nothing, and the
# old registration accepted it in silence: the contributing package's
# table then reads as complete while covering less than it claims,
# which is the worst thing a compatibility surface can do. It is not
# hypothetical. A drift-diffusion extension wrote `mixture()` where the
# vocabulary says `mixture`, and named `dec()` before that term was a
# feature at all, and lost two rows for a release without a word.
#
# Every test here restores both registries, because a registration
# lasts for the session and the rest of the suite reads the same table.

local_registries <- function(env = parent.frame()) {
  old_f <- frmtmb_compat_contrib$features
  old_r <- frmtmb_compat_contrib$rules
  old_a <- frmtmb_aterm_registry$reg
  withr::defer({
    frmtmb_compat_contrib$features <- old_f
    frmtmb_compat_contrib$rules <- old_r
    frmtmb_aterm_registry$reg <- old_a
  }, envir = env)
}

# a rule set of one row, built the way a contributor builds one
one_rule <- function(a, b, status = "works", note = "planted") {
  function() {
    bb <- compat_rule_builder()
    bb$r(a, b, status, note)
    bb$rules()
  }
}

test_that("a rule naming a feature outside the vocabulary is refused", {
  local_registries()
  err <- expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "mixture()")),
    class = "simpleError")
  msg <- conditionMessage(err)
  # the rule, the spelling, and the near miss the registrant meant
  expect_match(msg, "`wiener x mixture()`", fixed = TRUE)
  expect_match(msg, "'mixture()'", fixed = TRUE)
  expect_match(msg, "Did you mean 'mixture'?", fixed = TRUE)
  # and the two ways out
  expect_match(msg, "features =", fixed = TRUE)
  expect_match(msg, "expects =", fixed = TRUE)
})

test_that("a refused registration leaves the vocabulary as it found it", {
  local_registries()
  before <- frm_compat_features()
  expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "mixture()")))
  # the features of a rule set that never took would be a vocabulary
  # entry for a family the table says nothing about
  after <- frm_compat_features()
  expect_identical(after, before)
  expect_false("wiener" %in% after$name)
})

test_that("an unknown spelling with no near miss is refused without a guess", {
  local_registries()
  err <- expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "dec()")))
  # a confident wrong suggestion costs more than none
  expect_false(grepl("Did you mean", conditionMessage(err), fixed = TRUE))
  expect_match(conditionMessage(err), "'dec()'", fixed = TRUE)
})

test_that("a rule may name the features its own registration adds", {
  local_registries()
  expect_silent(
    frmtmb_register_compat(features = c(wiener = "family",
                                        "dec()" = "aterm"),
                           rules = one_rule("wiener", "dec()")))
  expect_equal(frm_compat("wiener", "dec()")$status, "works")
})

test_that("expects = admits another package's feature without adding it", {
  local_registries()
  expect_silent(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "hmm"),
                           expects = "hmm"))
  # the rule lies dormant: the feature is not in the table, because
  # nothing in this session implements it
  expect_false("hmm" %in% frm_compat_features()$name)
  expect_true("wiener" %in% frm_compat_features()$name)
})

test_that("a rule pattern names a real kind and a real group", {
  local_registries()
  expect_error(
    frmtmb_register_compat(rules = one_rule("gaussian", "kind:fmaily")),
    "no such kind of feature")
  expect_error(
    frmtmb_register_compat(rules = one_rule("gaussian", "group:cdfs")),
    "no such feature group")
  expect_silent(
    frmtmb_register_compat(rules = one_rule("kind:family", "group:cdf")))
})

test_that("a rule declares one of the five statuses", {
  local_registries()
  expect_error(
    frmtmb_register_compat(rules = one_rule("gaussian", "cens()", "refuse")),
    "declares the status 'refuse'", fixed = TRUE)
})

test_that("a rule naming one feature on both sides is refused", {
  local_registries()
  # the same failure as an unresolvable side, by another route: the
  # resolved table holds unordered pairs of DISTINCT features, so a
  # self-pair rule matches nothing however well spelled it is
  expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "wiener")),
    "names one feature on both sides")
  # a kind or a group on both sides covers distinct members and stands
  expect_silent(
    frmtmb_register_compat(rules = one_rule("kind:aterm", "kind:aterm")))
})

test_that("features = and rules = are checked for shape", {
  local_registries()
  expect_error(frmtmb_register_compat(features = c("family")),
               "DISPLAY name to its kind")
  expect_error(frmtmb_register_compat(features = c(wiener = "familly")),
               "no such kind")
  expect_error(frmtmb_register_compat(rules = data.frame(a = 1)),
               "FUNCTION of no arguments")
  expect_error(frmtmb_register_compat(rules = function() data.frame(a = 1)),
               "carrying feature_a")
  expect_error(frmtmb_register_compat(expects = 1),
               "as a character vector")
})

# ------------------------------------------- the addition-term seam

test_that("a registered addition term becomes a compat feature", {
  local_registries()
  expect_false("dec()" %in% frm_compat_features()$name)
  frmtmb_register_aterm("dec", arity = 1L)
  ft <- frm_compat_features()
  expect_equal(sum(ft$name == "dec()"), 1L)
  expect_equal(ft$kind[ft$name == "dec()"], "aterm")
  expect_equal(ft$key[ft$name == "dec()"], "dec")
  # so the table can be asked about it, rather than refusing a term
  # frm() accepts
  expect_equal(nrow(frm_compat("dec()", "gaussian")), 1L)
})

test_that("declaring a registered term again is a no-op, not a duplicate", {
  local_registries()
  # the order the drift-diffusion extension uses: the term first, then
  # the rules that name it, with the feature spelled out anyway
  frmtmb_register_aterm("dec", arity = 1L)
  expect_silent(
    frmtmb_register_compat(features = c(wiener = "family",
                                        "dec()" = "aterm"),
                           rules = one_rule("wiener", "dec()")))
  ft <- frm_compat_features()
  expect_equal(sum(ft$name == "dec()"), 1L)
  expect_equal(anyDuplicated(ft$name), 0L)
  expect_equal(frm_compat("wiener", "dec()")$status, "works")
})

test_that("one display name carries one kind", {
  local_registries()
  frmtmb_register_aterm("dec", arity = 1L)
  expect_error(
    frmtmb_register_compat(features = c("dec()" = "special")),
    "already holds that name under the kind 'aterm'", fixed = TRUE)
})

test_that("a name repeated in one features = call under two kinds is refused", {
  local_registries()
  # the cross-call clash one scope inward. Dropping the second entry
  # accepted a declaration that did nothing and said nothing, which is
  # the defect the whole seam exists to refuse
  expect_error(
    frmtmb_register_compat(features = c(zzq = "family", zzq = "aterm")),
    "the same call already gives that name the kind 'family'",
    fixed = TRUE)
  expect_false("zzq" %in% frm_compat_features()$name)
  # the same kind twice is the harmless half and stays silent: a rule
  # set may spell out a term frmtmb_register_aterm() already declared
  expect_silent(
    frmtmb_register_compat(features = c(zzq = "family", zzq = "family")))
  expect_equal(sum(frm_compat_features()$name == "zzq"), 1L)
})

test_that("an exact vocabulary entry beats an approximate match", {
  local_registries()
  # `se` used to be answered "Did you mean 'us'?": both halves of the
  # parenthesis test checked the same direction, so a name MISSING its
  # parentheses fell through to agrep(), which ranked us at one edit
  # above se() at two. An addition term pointed at an unstructured
  # covariance is the confident wrong suggestion the tight threshold
  # was chosen to avoid
  err <- expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "se")))
  expect_match(conditionMessage(err), "Did you mean 'se()'?", fixed = TRUE)
  expect_false(grepl("'us'", conditionMessage(err), fixed = TRUE))
  err2 <- expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "trials")))
  expect_match(conditionMessage(err2), "Did you mean 'trials()'?",
               fixed = TRUE)
})

test_that("a rules = builder that throws says which argument failed", {
  local_registries()
  # the builder is called HERE now, so an error in it takes down the
  # registrant's .onLoad(). Raw, it arrived out of loadNamespace() with
  # nothing naming the argument that called it
  err <- expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = function() stop("boom")))
  msg <- conditionMessage(err)
  expect_match(msg, "frmtmb_register_compat()", fixed = TRUE)
  expect_match(msg, "rules = builder", fixed = TRUE)
  expect_match(msg, "boom", fixed = TRUE)
  expect_false("wiener" %in% frm_compat_features()$name)
})

test_that("the failure clause names the package whose load it was", {
  # registration runs from a contributor's .onLoad(), where the user has
  # no call stack to read
  expect_identical(compat_registrant(asNamespace("stats")),
                   " that stats gave it")
  expect_identical(compat_registrant(globalenv()), "")
  expect_identical(compat_registrant(asNamespace("frmtmb")), "")
})

test_that("expects = refuses a name the vocabulary already holds", {
  local_registries()
  # expects = exempts a name from the resolvability check, so a name the
  # session already has exempts nothing
  expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "cens()"),
                           expects = "cens()"),
    "exempts nothing")
  expect_false("wiener" %in% frm_compat_features()$name)
  # including one this very call supplies
  expect_error(
    frmtmb_register_compat(features = c(wiener = "family"),
                           rules = one_rule("wiener", "cens()"),
                           expects = "wiener"),
    "exempts nothing")
})

test_that("a declared but unresolved rule side is reported, not dropped", {
  local_registries()
  expect_null(attr(frm_compat("gaussian", "cens()"), "unresolved"))
  frmtmb_register_compat(features = c(wiener = "family"),
                         rules = one_rule("wiener", "hmm"),
                         expects = "hmm")
  out <- frm_compat("wiener")
  # it cannot honestly be a row: nothing here implements hmm, so there
  # is no pair. Saying nothing is what expects = re-opened, so the
  # resolver says it at the one moment the session is known
  expect_false(any(out$feature_b == "hmm"))
  expect_identical(attr(out, "unresolved"), "hmm")
  expect_s3_class(out, "frmtmb_compat")
  expect_output(print(out), "Unresolved rule sides")
  expect_output(print(out), "hmm")
  # and goes quiet the moment something supplies it
  frmtmb_register_compat(features = c(hmm = "structure"))
  out2 <- frm_compat("wiener")
  expect_null(attr(out2, "unresolved"))
  expect_false(inherits(out2, "frmtmb_compat"))
  expect_equal(frm_compat("wiener", "hmm")$status, "works")
})

test_that("a refused frmtmb_register_aterm() registers nothing at all", {
  local_registries()
  before_ft <- frm_compat_features()
  before_reg <- names(frmtmb_aterm_registry$reg)
  # thirteen vocabulary entries are written with parentheses and are not
  # of kind aterm. The parser used to be committed first, so the call
  # threw AND registered: frm() then took `y | s(col) ~ x` while
  # frm_compat() described s() as a smooth, which is the split state
  # this seam exists to make unreachable
  for (nm in c("s", "ar", "mm")) {
    expect_error(frmtmb_register_aterm(nm),
                 "already means something else in a formula")
    expect_false(nm %in% names(frmtmb_aterm_registry$reg))
  }
  expect_identical(names(frmtmb_aterm_registry$reg), before_reg)
  expect_identical(frm_compat_features(), before_ft)
  # and the parser is where the consequence would have shown
  expect_null(registered_aterm_of("s"))
})

test_that("the refusal reads in the term registrant's own terms", {
  local_registries()
  msg <- conditionMessage(expect_error(frmtmb_register_aterm("s")))
  expect_match(msg, "frmtmb_register_aterm()", fixed = TRUE)
  expect_match(msg, "s()", fixed = TRUE)
  expect_match(msg, "a 'special' feature", fixed = TRUE)
  expect_match(msg, "nothing has been registered", fixed = TRUE)
})

test_that("a kind beginning with a vowel takes the right article", {
  local_registries()
  frmtmb_register_aterm("dec", arity = 1L)
  expect_match(
    conditionMessage(expect_error(
      frmtmb_register_compat(features = c("dec()" = "special")))),
    "as a 'special' feature", fixed = TRUE)
  expect_match(
    conditionMessage(expect_error(
      frmtmb_register_compat(features = c(mixture = "aterm")))),
    "as an 'aterm' feature", fixed = TRUE)
})

