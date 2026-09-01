# Tests for the feature compatibility registry.
#
# Two jobs. First, keep the registry internally sound and anchored to
# the package's real vocabulary, so a renamed or removed feature cannot
# leave a stale row behind. Second, check a sample of "refused" claims
# against what frm() actually does, so the registry cannot drift into
# describing refusals that no longer happen.
#
# Only refusals whose guards exist today are exercised. "broken" rows
# are deliberately not run: they record defects, and asserting the
# defect would freeze it in place.

# ------------------------------------------------------------ invariants

test_that("registry statuses are drawn from the declared set", {
  rules <- frm_compat_rules()
  expect_true(all(rules$status %in% frmtmb_compat_statuses))
  expect_true(all(frm_compat()$status %in% frmtmb_compat_statuses))
})

test_that("feature names are unique and every kind is populated", {
  ft <- frm_compat_features()
  expect_equal(anyDuplicated(ft$name), 0L)
  expect_true(all(nzchar(ft$name)), all(nzchar(ft$key)))
  expect_setequal(
    unique(ft$kind),
    c("family", "covstruct", "aterm", "special", "mode", "structure",
      "method", "grammar"))
})

test_that("no rule is stated twice for the same unordered pattern pair", {
  rules <- frm_compat_rules()
  key <- vapply(seq_len(nrow(rules)), function(i) {
    paste(sort(c(rules$feature_a[i], rules$feature_b[i])),
          collapse = " <-> ")
  }, "")
  dup <- key[duplicated(key)]
  # A duplicate silently shadows the earlier rule, which is exactly the
  # kind of invisible override the registry exists to prevent.
  expect_equal(dup, character(0))
})

test_that("rule patterns name real kinds, groups, and features", {
  ft <- frm_compat_features()
  rules <- frm_compat_rules()
  pats <- unique(c(rules$feature_a, rules$feature_b))
  for (p in pats) {
    if (p == "*") next
    if (startsWith(p, "kind:")) {
      expect_true(substring(p, 6L) %in% ft$kind, info = p)
    } else if (startsWith(p, "group:")) {
      expect_true(substring(p, 7L) %in% names(frmtmb_compat_groups_lst),
                  info = p)
    } else {
      expect_true(p %in% ft$name, info = p)
    }
  }
})

test_that("every group member is a declared feature", {
  ft <- frm_compat_features()
  for (g in names(frmtmb_compat_groups_lst)) {
    expect_true(all(frmtmb_compat_groups_lst[[g]] %in% ft$name),
                info = g)
  }
})

test_that("refused and broken rows carry an explanation", {
  res <- frm_compat(status = c("refused", "broken"))
  expect_true(all(nzchar(res$note)))
})

test_that("every resolved pair gets a status", {
  res <- frm_compat()
  expect_false(any(is.na(res$status)))
  # family x family pairs are meaningless; a model carries one family
  expect_false(any(res$kind_a == "family" & res$kind_b == "family"))
})

# ------------------------------------- the registry names real package parts

test_that("declared families exist in the family registry", {
  ft <- frm_compat_features()
  fams <- ft$key[ft$kind == "family"]
  expect_true(all(fams %in% names(family_registry)))
})

test_that("declared covariance structures exist in the covstruct registry", {
  ft <- frm_compat_features()
  covs <- ft$key[ft$kind == "covstruct"]
  expect_true(all(covs %in% names(covstruct_registry)))
})

test_that("declared addition terms are the ones the parser accepts", {
  ft <- frm_compat_features()
  ats <- ft$key[ft$kind == "aterm"]
  d <- data.frame(y = rnorm(10), x = rnorm(10))
  # the parser lists its whole supported set when it rejects one
  msg <- tryCatch(
    frm(y | not_an_aterm(x) ~ x, data = d, family = gaussian(),
        dry_run = "spec"),
    error = conditionMessage)
  for (a in ats) {
    expect_true(grepl(paste0(a, "()"), msg, fixed = TRUE), info = a)
  }
})

test_that("declared specials still parse", {
  ft <- frm_compat_features()
  expect_setequal(ft$key[ft$kind == "special"],
                  c("s", "t2", "mo", "mi", "gp", "cs"))
  set.seed(1)
  n <- 60
  d <- data.frame(
    y = rnorm(n), x = rnorm(n), x2 = rnorm(n),
    o = factor(sample(1:4, n, TRUE), ordered = TRUE),
    g = factor(rep(1:6, each = 10)))
  # dry_run stops after the model frame, so this stays cheap
  expect_s3_class(frm(y ~ s(x, k = 5), data = d, family = gaussian(),
                      dry_run = "frame"), "frmtmb_frame")
  expect_s3_class(frm(y ~ t2(x, x2), data = d,
                      family = gaussian(), dry_run = "frame"),
                  "frmtmb_frame")
  expect_s3_class(frm(y ~ mo(o), data = d, family = gaussian(),
                      dry_run = "frame"), "frmtmb_frame")
  expect_s3_class(frm(y ~ gp(x, k = 5), data = d, family = gaussian(),
                      dry_run = "frame"), "frmtmb_frame")
  expect_s3_class(frm(o ~ cs(x), data = d, family = sratio(),
                      dry_run = "frame"), "frmtmb_frame")
})

test_that("declared modes are arguments of frm() or frmtmb_control()", {
  ft <- frm_compat_features()
  modes <- ft$key[ft$kind == "mode"]
  known <- c(names(formals(frm)), names(formals(frmtmb_control)))
  # "bounds" is the registry's name for the lower/upper argument pair
  known <- c(known, if (all(c("lower", "upper") %in% known)) "bounds")
  expect_true(all(modes %in% known))
})

test_that("declared post-fit methods are callable", {
  ft <- frm_compat_features()
  meths <- ft$key[ft$kind == "method"]
  # Several rows name an argument setting of a generic, or a method
  # registered on a suggested package's generic, rather than a
  # function of their own.
  as_fn <- c(residuals_osa = "residuals.frmtmb_fit",
             confint_profile = "confint.frmtmb_fit",
             hypothesis_profile = "hypothesis",
             emmeans = "emm_basis.frmtmb_fit",
             fitted = "fitted.frmtmb_fit",
             predict = "predict.frmtmb_fit",
             simulate = "simulate.frmtmb_fit",
             residuals = "residuals.frmtmb_fit")
  for (m in meths) {
    fn <- if (m %in% names(as_fn)) as_fn[[m]] else m
    expect_true(exists(fn, mode = "function"), info = m)
  }
})

# ------------------------------------------------------------ query surface

test_that("frm_compat() slices the table three ways", {
  all_pairs <- frm_compat()
  expect_gt(nrow(all_pairs), 1000)
  expect_named(all_pairs,
               c("feature_a", "kind_a", "feature_b", "kind_b",
                 "status", "note"))

  one <- frm_compat("trunc()")
  expect_true(all(one$feature_a == "trunc()"))
  expect_lt(nrow(one), nrow(all_pairs))

  pair <- frm_compat("rescor", "cens()")
  expect_equal(nrow(pair), 1L)
  expect_equal(pair$status, "refused")
  # unordered: the answer does not depend on the argument order
  expect_equal(frm_compat("cens()", "rescor")$status, "refused")
})

test_that("frm_compat() filters by status and rejects unknown input", {
  broken <- frm_compat(status = "broken")
  expect_true(all(broken$status == "broken"))
  expect_gt(nrow(broken), 0L)
  expect_error(frm_compat("no_such_feature"), "Unknown feature")
  expect_error(frm_compat(status = "maybe"), "Unknown status")
})

test_that("a more specific rule beats a broader one", {
  # kind:family x cens() refuses, group:cdf_continuous grants, and the
  # explicit poisson row refuses again
  expect_equal(frm_compat("cens()", "gaussian")$status, "works")
  expect_equal(frm_compat("cens()", "poisson")$status, "refused")
  expect_equal(frm_compat("cens()", "beta")$status, "refused")
  expect_equal(frm_compat("trunc()", "poisson")$status, "conditional")
  expect_equal(frm_compat("trunc()", "weibull")$status, "works")
})

# ----------------------------------------- declared refusals really refuse

test_that("declared refusals at the frame stage really refuse", {
  set.seed(2)
  n <- 60
  d <- data.frame(
    y = rnorm(n, 5), y2 = rnorm(n), x = rnorm(n),
    s = runif(n, .2, .5), w = runif(n, .5, 2),
    cc = rep(0L, n),
    o = factor(sample(1:3, n, TRUE), ordered = TRUE),
    g = factor(rep(1:6, each = 10)))
  fr <- function(...) frm(..., data = d, dry_run = "frame")

  refuse <- function(a, b, expr) {
    expect_equal(frm_compat(a, b)$status, "refused")
    expect_error(expr)
  }
  refuse("rescor", "cens()",
         fr(bf(y | cens(cc) ~ x) + bf(y2 ~ x) + set_rescor(TRUE),
            family = gaussian()))
  refuse("rescor", "se()",
         fr(bf(y | se(s) ~ x) + bf(y2 ~ x) + set_rescor(TRUE),
            family = gaussian()))
  refuse("rescor", "weights()",
         fr(bf(y | weights(w) ~ x) + bf(y2 ~ x) + set_rescor(TRUE),
            family = gaussian()))
  refuse("rescor", "poisson",
         fr(bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE),
            family = poisson()))
  refuse("se()", "poisson",
         fr(y | se(s) ~ x, family = poisson()))
  refuse("cens()", "binomial",
         fr(y | cens(cc) ~ x, family = bernoulli()))
  refuse("cs_pred()", "cumulative",
         fr(o ~ cs(x), family = cumulative()))
  refuse("|ID|", "ar1",
         fr(bf(y ~ ar1(0 + o | q | g)) + bf(y2 ~ x) + set_rescor(FALSE),
            family = gaussian()))
})

test_that("declared refusals at the fit stage really refuse", {
  set.seed(3)
  n <- 60
  d <- data.frame(y = rnorm(n, 5), x = rnorm(n),
                  g = factor(rep(1:6, each = 10)),
                  o = factor(rep(1:10, 6), ordered = TRUE))

  expect_equal(frm_compat("REML", "quadrature")$status, "refused")
  expect_error(frm(y ~ x + (1 | g), data = d, family = gaussian(),
                   REML = TRUE, quadrature = TRUE),
               "quadrature")

  expect_equal(frm_compat("REML", "profile")$status, "refused")
  expect_error(frm(y ~ x + (1 | g), data = d, family = gaussian(),
                   REML = TRUE,
                   control = frmtmb_control(profile = TRUE)),
               "profile")

  expect_equal(frm_compat("quadrature", "ar1")$status, "refused")
  expect_error(frm(y ~ ar1(0 + o | g), data = d, family = gaussian(),
                   quadrature = TRUE),
               "scalar random")

  expect_equal(frm_compat("profile", "bounds")$status, "refused")
  expect_error(frm(y ~ x, data = d, family = gaussian(),
                   lower = c("x" = -10),
                   control = frmtmb_control(profile = TRUE)),
               "bounds")
})
