# The pre-flight: frm_sample() reads the compatibility registry for
# rows that refuse it, and stops before anything is taped.
#
# A refusal declared in the registry and enforced nowhere is a note, not
# a refusal. This file is the other half of `?frmtmb_register_compat`'s
# promise that a `refused` row is consulted.

skip_on_cran()

pf_data <- function(seed = 4, n = 40L, ng = 5L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x, 1)
  dd
}

## ---- what the pre-flight can see -------------------------------------

test_that("formula_calls() finds a call wherever a model writes one", {
  f <- function(...) frmtmb.sample:::formula_calls(...)$other
  expect_true("s" %in% f(bf(y ~ s(x) + (1 | g))))
  # a nonlinear body, which is where frm_ode() sits
  expect_true("frm_ode" %in%
                f(bf(y ~ frm_ode(dyn, init = list(d), times = t),
                     a ~ 1, nl = TRUE)))
  # an nlf() body reached through nlforms
  expect_true("frm_ode" %in%
                f(bf(y ~ a * b, a ~ 1, nl = TRUE) +
                    nlf(b ~ frm_ode(dyn, init = list(d), times = t))))
  # every response of a multivariate model
  expect_true("s" %in% f(mvbf(bf(y1 ~ x), bf(y2 ~ s(x)))))
  # and a model that writes none is not credited with one
  expect_false("frm_ode" %in% f(bf(y ~ x + (1 | g))))
})

test_that("formula_calls() puts an addition term in aterm position and
           the same name elsewhere in `other`", {
  f <- frmtmb.sample:::formula_calls
  # the split that separates the se() ADDITION TERM from a user's own
  # se() in a nonlinear body, which frmtmb fits and must not refuse
  a <- f(bf(y | se(v) ~ x))
  expect_true("se" %in% a$aterm)
  expect_false("se" %in% a$other)
  b <- f(bf(y ~ a * se(x), a ~ 1, nl = TRUE))
  expect_false("se" %in% b$aterm)
  expect_true("se" %in% b$other)
  # several addition terms at once
  cc <- f(bf(y | cens(c) + weights(w) ~ x))
  expect_true(all(c("cens", "weights") %in% cc$aterm))
  # the RESPONSE expression is not aterm position, so cbind() does not
  # become one
  d <- f(bf(cbind(y1, y2) | trials(n) ~ x))
  expect_true("trials" %in% d$aterm)
  expect_false("cbind" %in% d$aterm)
  expect_true("cbind" %in% d$other)
  # a `|` inside a dpar formula is a random-effect bar, not an aterm
  e <- f(bf(y ~ x, sigma ~ x + (1 | g)))
  expect_length(e$aterm, 0L)
  # and a predictor special is never aterm position
  g <- f(bf(y ~ s(x) + (1 | id)))
  expect_length(g$aterm, 0L)
})

test_that("model_family_names() reads both routes", {
  g <- frmtmb.sample:::model_family_names
  expect_equal(g(bf(y ~ x) + poisson()), "poisson")
  expect_true("poisson" %in% g(bf(y ~ x), family = poisson()))
  expect_setequal(g(mvbf(bf(y1 ~ x) + poisson(), bf(y2 ~ x) + gaussian())),
                  c("poisson", "gaussian"))
  fit <- frm(bf(y ~ x) + gaussian(), data = pf_data())
  expect_equal(g(fit), "gaussian")
})

## ---- the three branches of the decision ------------------------------

# The registry is global and appends for the life of the session, so
# these drive sample_preflight() with rows supplied directly rather than
# registering a synthetic feature that every later file would inherit.
# The registry READ is covered end to end by the frmtmb.ode block below.

pf_rows <- function(feature, kind, note = "Refused for the test.") {
  data.frame(feature_a = "frm_sample", kind_a = "method",
             feature_b = feature, kind_b = kind,
             status = "refused", note = note,
             stringsAsFactors = FALSE)
}

test_that("a refused call-shaped feature the model writes stops the
           call, and one it does not write lets it through", {
  pf <- frmtmb.sample:::sample_preflight
  rows <- pf_rows("frm_ode()", "special", "Broken at the boundary.")
  form <- bf(y ~ frm_ode(dyn, init = list(d), times = t), a ~ 1,
             nl = TRUE)
  expect_error(pf(form, rows = rows),
               "frm_ode[(][)].*Broken at the boundary")
  expect_silent(pf(bf(y ~ x + (1 | g)), rows = rows))
})

test_that("a refused family is matched on the family name", {
  pf <- frmtmb.sample:::sample_preflight
  rows <- pf_rows("poisson", "family", "No.")
  expect_error(pf(bf(y ~ x) + poisson(), rows = rows), "poisson")
  expect_error(pf(bf(y ~ x), family = poisson(), rows = rows), "poisson")
  expect_silent(pf(bf(y ~ x) + gaussian(), rows = rows))
})

test_that("a refused row the pre-flight cannot match warns rather than
           passing in silence", {
  pf <- frmtmb.sample:::sample_preflight
  # a fitting mode is neither a call in a formula nor a family name
  rows <- pf_rows("REML", "mode", "Not supported.")
  expect_warning(pf(bf(y ~ x) + gaussian(), rows = rows),
                 "cannot tell whether this model uses it")
  # the reason is part of the warning, because "cannot tell" without a
  # reason gives a registrant nothing to act on
  expect_warning(pf(bf(y ~ x) + gaussian(), rows = rows),
                 "without parentheses")
})

## ---- deciding by position, and refusing to decide by name -----------
#
# What decides whether a model uses a feature is WHERE a call sits in
# the formula. One position is available without reading anything from
# core, the addition-term slot, and these blocks pin what it settles and
# what is left over. See dev/stanctl-findings.md for the seam that would
# settle the rest.

test_that("an addition term is decided by position, so the same call
           in a nonlinear body does not refuse", {
  pf <- frmtmb.sample:::sample_preflight
  rows <- pf_rows("se()", "aterm", "No.")
  # the ADDITION TERM: refused
  expect_error(pf(bf(y | se(v) ~ x), rows = rows), "se[(][)]")
  # the user's own helper of the same name in a nonlinear body. frmtmb
  # FITS this model, so refusing it would be a check firing on a
  # correct model, which the lane rules rank below no check at all.
  expect_silent(pf(bf(y ~ a * se(x), a ~ 1, nl = TRUE), rows = rows))
  # position, not the name, is doing the work: neither a uniqueness nor
  # a base R test would separate these two, and `trunc()` is an aterm
  # that IS a base R function
  rows2 <- pf_rows("trunc()", "aterm", "No.")
  expect_error(pf(bf(y | trunc(lb = 0) ~ x), rows = rows2), "trunc")
  expect_silent(pf(bf(y ~ a * trunc(x), a ~ 1, nl = TRUE), rows = rows2))
})

test_that("the mi()/mi_pred() split is settled in both directions", {
  pf <- frmtmb.sample:::sample_preflight
  ft <- frm_compat_features()
  # core gives each of three pairs ONE parser key while splitting them
  # by display name (R/compat.R, "The special vocabulary collides with
  # the covariance vocabulary"). This is the drift guard for that.
  expect_equal(sort(unique(ft$key[duplicated(ft$key)])),
               c("cs", "gp", "mi"))
  # the aterm row fires on the aterm and not on the predictor special
  a <- pf_rows("mi()", "aterm", "No.")
  expect_error(pf(bf(y | mi() ~ x), rows = a), "mi[(][)]")
  expect_silent(pf(bf(y ~ mi(x)), rows = a))
  # and the special row decides neither, because no formula writes the
  # call `mi_pred`. It says so rather than passing quietly; the block
  # below on display names a formula does not write is where that is
  # pinned, and the filed seam is what would decide it.
  b <- pf_rows("mi_pred()", "special", "No.")
  expect_warning(pf(bf(y | mi() ~ x), rows = b), "no formula writes")
  expect_warning(pf(bf(y ~ mi(x)), rows = b), "no formula writes")
})

test_that("a name-route feature that is also a base R function is not
           decided by guessing", {
  pf <- frmtmb.sample:::sample_preflight
  brf <- frmtmb.sample:::base_r_function
  expect_true(brf("trunc"))
  expect_true(brf("weights"))
  expect_true(brf("ar"))
  expect_false(brf("frm_ode"))
  # all fourteen base-priority packages, not thirteen
  expect_true(brf("tclvalue") || !requireNamespace("tcltk",
                                                   quietly = TRUE))
  # `ar()` is an autocor, so it has no position of its own and falls to
  # the name route, where the base R test stops it
  rows <- pf_rows("ar()", "autocor", "No.")
  expect_warning(pf(bf(y ~ a * ar(x), a ~ 1, nl = TRUE), rows = rows),
                 "base R function")
  # a bare display name is the vocabulary saying the feature is not
  # written as a call, so it takes no route at all
  rows2 <- pf_rows("exp", "covstruct", "No.")
  expect_warning(pf(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE),
                    rows = rows2),
                 "without parentheses")
})

test_that("the premises the un-revocable match rests on", {
  # Matching the DISPLAY name rather than the parser key is what makes a
  # registered refusal un-revocable, and the closure is a proof rather
  # than a measurement. These are its premises, asserted so a change to
  # core that breaks one fails here rather than in the field.
  ft <- frm_compat_features()
  # 1. core refuses a second registration of a display name, so display
  #    names are unique. The KIND is as load-bearing as the name: the
  #    route is chosen from both, and core's refusal covers both.
  expect_equal(sum(duplicated(ft$name)), 0L)
  expect_error(
    frmtmb_register_compat(features = c("s()" = "grammar")),
    "One display name carries one kind")
  # 2. dropping a fixed "()" suffix is injective over the names that
  #    carry it, so no two call-shaped rows want the same call
  callable <- ft$name[grepl("[(][)]$", ft$name)]
  expect_gt(length(callable), 0L)
  expect_equal(sum(duplicated(sub("[(][)]$", "", callable))), 0L)
  # 3. a bare name cannot reach the call route. This is a property of
  #    the CODE, not of the vocabulary: a formula DOES write `ar1(...)`
  #    and `gr(...)`, and what stops those rows is that
  #    preflight_route() returns the without-parentheses sentence
  #    before a bare name gets that far.
  ar1_calls <- frmtmb.sample:::formula_calls(
    bf(y ~ x + (1 | ar1(t, g))))
  expect_true("ar1" %in% ar1_calls$other)
  expect_match(frmtmb.sample:::preflight_route("ar1", "covstruct", ft),
               "without parentheses")
})

test_that("the bare-name collision that used to revoke the ODE refusal
           no longer does", {
  skip_if_not_installed("frmtmb.ode")
  loadNamespace("frmtmb.ode")
  route <- frmtmb.sample:::preflight_route
  form <- bf(y ~ frm_ode(d, init = list(1), times = t), a ~ 1,
             nl = TRUE)
  rows <- pf_rows("frm_ode()", "special", "No.")
  expect_equal(route("frm_ode()", "special", frm_compat_features()),
               "call")
  expect_error(frmtmb.sample:::sample_preflight(form, rows = rows),
               "frm_ode[(][)]")

  # The attack, performed for real rather than against a hand-built
  # table: core ACCEPTS a bare `frm_ode` display name, and under the old
  # parser-key match that shared key downgraded the refusal to a warning
  # and the model sampled. The registration is permanent for this R
  # process and there is no way to withdraw one, so it comes last in
  # this block; the blocks after it are unaffected because the bare row
  # carries no rules and so changes no pair's status.
  expect_silent(
    frmtmb_register_compat(features = c("frm_ode" = "structure")))
  ft2 <- frm_compat_features()
  expect_true(all(c("frm_ode()", "frm_ode") %in% ft2$name))
  expect_equal(length(unique(ft2$key[ft2$name %in%
                                       c("frm_ode()", "frm_ode")])), 1L)
  # the route does not move and the refusal still fires
  expect_equal(route("frm_ode()", "special", ft2), "call")
  expect_error(frmtmb.sample:::sample_preflight(form, rows = rows),
               "frm_ode[(][)]")
})

test_that("a display name a formula does not write warns instead of
           missing in silence", {
  # The price of matching the display name is the three `_pred`
  # specials: `mi_pred()` is spelled `mi` by a formula, so a name match
  # on `mi_pred` can never fire. Missing in silence there would be the
  # failing-open shape this file exists to close, and the registry
  # itself says which rows are affected, since its key IS the parser
  # name.
  pf <- frmtmb.sample:::sample_preflight
  route <- frmtmb.sample:::preflight_route
  ft <- frm_compat_features()
  differ <- ft$name[grepl("[(][)]$", ft$name) &
                      sub("[(][)]$", "", ft$name) != ft$key]
  expect_setequal(differ, c("mi_pred()", "gp_pred()", "cs_pred()"))
  for (nm in differ) {
    expect_match(route(nm, "special", ft), "no formula writes the call")
  }
  rows <- pf_rows("mi_pred()", "special", "No.")
  expect_warning(pf(bf(y ~ mi(x)), rows = rows),
                 "no formula writes the call")
  expect_warning(pf(bf(y | mi() ~ x), rows = rows),
                 "no formula writes the call")
  # and every other call-shaped row still decides, so the guard is not
  # a blanket
  expect_equal(route("s()", "special", ft), "call")
})

## ---- end to end, on the row frmtmb.ode registers ---------------------

test_that("frmtmb.ode registers a refused row and frm_sample() reads
           it before assembling anything", {
  skip_if_not_installed("frmtmb.ode")
  skip_sampler()
  loadNamespace("frmtmb.ode")
  rows <- frm_compat("frm_sample", status = "refused")
  expect_true("frm_ode()" %in% rows$feature_b)

  form <- bf(conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                            parms = list(exp(lka)), group = id),
             lka ~ 1, nl = TRUE) + gaussian()
  # `data` is not supplied. The formula route refuses a missing `data`
  # the moment it tries to assemble, so getting the ODE refusal instead
  # PROVES the pre-flight ran before assembly, and therefore before any
  # tape was built. No timing is involved.
  err <- expect_error(frm_sample(form, chains = 1, iter = 20))
  expect_match(conditionMessage(err), "frm_ode[(][)]")
  expect_false(grepl("needs data", conditionMessage(err)))
})

test_that("an ODE fit is refused on the fit route too, with no chain
           started", {
  skip_if_not_installed("frmtmb.ode")
  skip_if_not_installed("RTMBode")
  skip_sampler()
  loadNamespace("frmtmb.ode")
  pk_dyn <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(-p[1L] * y[1L], p[1L] * y[1L] - p[2L] * y[2L]))
  }
  dd <- data.frame(id = factor(rep(1:3, each = 4L)),
                   time = rep(c(0.5, 1, 2, 4), 3L), dose = 100,
                   conc = rep(c(3, 4, 3, 2), 3L))
  form <- bf(conc ~ frmtmb.ode::frm_ode(pk_dyn, init = list(dose, 0),
                                        times = time,
                                        parms = list(exp(lka), exp(lke)),
                                        group = id, output = 2L),
             lka ~ 1, lke ~ 1, nl = TRUE) + gaussian()
  fit <- frm(form, data = dd, dry_run = "objective",
             start = list(beta = c(0, log(0.25))))
  expect_error(frm_sample(fit, chains = 1, iter = 20, refresh = 0),
               "frm_ode[(][)]")
  # the qualified spelling is what the fit above uses, so the walk has
  # to see through `pkg::fn` as well as a bare call
  expect_true("frm_ode" %in%
                frmtmb.sample:::formula_calls(fit[["bform"]])$other)
  # as_tmbstan() is the OTHER door, and the one where the failure is
  # worse: it hands back an empty stanfit with no error at all, and the
  # abort behind it leaves rstan's autodiff arena unbalanced for the
  # rest of the session, so every later chain in that session fails too
  expect_error(as_tmbstan(fit, chains = 1, iter = 20, refresh = 0),
               "frm_ode[(][)]")
  expect_error(as_tmbstan(fit, chains = 1, iter = 20, refresh = 0),
               "as_tmbstan[(][)]")
})

test_that("the pre-flight leaves the conditional rows conditional", {
  skip_if_not_installed("frmtmb.ode")
  loadNamespace("frmtmb.ode")
  tb <- frm_compat("frm_sample")
  # the one refusal is the ODE row, and nothing else became one
  expect_equal(sort(tb$feature_b[tb$status == "refused"]), "frm_ode()")
})

test_that("the eam and latent rows stay conditional", {
  # This is the plan's second check for item 1.2, and it needs the
  # packages that OWN those features: wiener is frmtmb.eam's and hmm
  # and lca are frmtmb.latent's, and the suite runs one file per
  # process with neither attached. Written as a `next` over names the
  # vocabulary might not carry, this block ran zero assertions and
  # passed, which is the failing-open shape the project keeps closing.
  # It now loads what it needs, skips with a reason when it cannot, and
  # counts what it actually checked.
  owners <- c(wiener = "frmtmb.eam", hmm = "frmtmb.latent",
              lca = "frmtmb.latent")
  for (pkg in unique(owners)) {
    skip_if_not_installed(pkg)
    loadNamespace(pkg)
  }
  skip_if_not_installed("frmtmb.ode")
  loadNamespace("frmtmb.ode")
  tb <- frm_compat("frm_sample")
  vocab <- frm_compat_features()$name
  checked <- 0L
  for (nm in names(owners)) {
    expect_true(nm %in% vocab)
    st <- tb$status[tb$feature_b == nm]
    expect_length(st, 1L)
    expect_false(identical(st, "refused"))
    checked <- checked + 1L
  }
  # the guard on the guard: a future version of this block that visits
  # nothing fails here instead of passing
  expect_equal(checked, length(owners))
})
