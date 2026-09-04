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
    c("family", "covstruct", "aterm", "autocor", "special", "mode",
      "structure", "method", "grammar"))
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

test_that("no rule names the same feature on both sides", {
  rules <- frm_compat_rules()
  # The resolved table holds unordered pairs of DISTINCT features, so a
  # rule naming one feature twice can never match. Two rules may share
  # a PATTERN on both sides ("kind:covstruct" x "kind:covstruct"),
  # which covers distinct members of that kind.
  self <- rules$feature_a == rules$feature_b &
    compat_spec(rules$feature_a) == 3L
  expect_equal(rules$feature_a[self], character(0))
})

test_that("every rule wins at least one pair", {
  pairs <- frmtmb_compat_pairs_tbl()
  rules <- frm_compat_rules()
  win <- compat_resolve(pairs, rules)$win
  dead <- setdiff(seq_len(nrow(rules)), unique(win))
  # An unreachable rule is a claim nobody can read: it looks like a
  # declaration and behaves like a comment.
  expect_equal(
    sprintf("%s x %s", rules$feature_a[dead], rules$feature_b[dead]),
    character(0))
})

test_that("no pair is decided by file order alone", {
  # Two rules of the same specificity signature that disagree about a
  # pair leave the answer to whichever was typed last. Every such tie
  # must be a declared override.
  bad <- frmtmb_compat_validate()
  expect_equal(sprintf("%s x %s", bad$feature_a, bad$feature_b),
               character(0))
})

test_that("precedence compares the two sides as a sorted pair", {
  # (3,1) beats (2,2) beats (2,1) beats (1,1); summing the sides would
  # tie the first two at 4 and the last two at 3.
  expect_gt(compat_sig_rank("cens()", "kind:family"),
            compat_sig_rank("group:cdf", "group:discrete"))
  expect_gt(compat_sig_rank("group:cdf", "group:discrete"),
            compat_sig_rank("group:cdf", "kind:family"))
  expect_gt(compat_sig_rank("group:cdf", "kind:family"),
            compat_sig_rank("kind:family", "kind:aterm"))
  # a rule more specific on one side and no less specific on the other
  # always wins, whichever side that is
  expect_gt(compat_sig_rank("cens()", "*"),
            compat_sig_rank("group:cdf", "*"))
  expect_gt(compat_sig_rank("*", "cens()"),
            compat_sig_rank("*", "group:cdf"))
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

test_that("declared autocorrelation terms parse", {
  ft <- frm_compat_features()
  expect_setequal(ft$key[ft$kind == "autocor"],
                  c("ar", "ma", "arma", "cosy", "unstr"))
  set.seed(1)
  d <- expand.grid(week = 1:4, subj = factor(1:8))
  d$y <- rnorm(32)
  d$x <- rnorm(32)
  for (tm in c("ar(week, subj, cov = TRUE)", "ma(week, subj, cov = TRUE)",
               "arma(week, subj, cov = TRUE)", "cosy(week, subj)",
               "unstr(week, subj)")) {
    fr <- frm(stats::as.formula(paste("y ~ x +", tm)), data = d,
              family = gaussian(), dry_run = "frame")
    expect_s3_class(fr, "frmtmb_frame")
    expect_length(fr$autocor, 1L)
  }
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

test_that("frm_compat() crosses vector arguments", {
  # a vector used to be recycled against the whole table, which gave a
  # silently wrong slice rather than the cross the caller asked for
  many <- frm_compat(c("cens()", "trunc()"), c("gaussian", "poisson"))
  expect_equal(nrow(many), 4L)
  expect_setequal(many$feature_a, c("cens()", "trunc()"))
  expect_setequal(many$feature_b, c("gaussian", "poisson"))
  expect_equal(many$status[many$feature_a == "cens()" &
                             many$feature_b == "poisson"], "refused")

  one_side <- frm_compat(c("cens()", "trunc()"))
  expect_true(all(one_side$feature_a %in% c("cens()", "trunc()")))
  expect_equal(nrow(one_side),
               nrow(frm_compat("cens()")) + nrow(frm_compat("trunc()")) - 1L)
})

test_that("frm_compat() rejects an empty feature argument", {
  # returning an empty table reads exactly like "this feature takes
  # part in no pairs", which is never true
  expect_error(frm_compat(character(0)), "empty")
  expect_error(frm_compat("cens()", character(0)), "empty")
  expect_error(frm_compat(42), "character")
})

test_that("frm_compat() filters by status and rejects unknown input", {
  # the broken bucket can be (and currently is) empty: v0.24 cleared
  # every known-broken pair; the filter must still work on it
  broken <- frm_compat(status = "broken")
  expect_true(all(broken$status == "broken"))
  refused <- frm_compat(status = "refused")
  expect_true(all(refused$status == "refused"))
  expect_gt(nrow(refused), 0L)
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

test_that("a structure's own condition survives the family", {
  # A <covstruct> x kind:family rule granting "works" used to outrank
  # each structure's own condition, so 429 pairs claimed unconditional
  # support for a spatial or positional block. exp() over a plain
  # factor still builds, with the level order standing in for the
  # coordinates - which is the condition, not a refusal.
  for (cs in c("exp", "gau", "mat", "ou", "ar1", "rr", "gr_cov")) {
    expect_equal(frm_compat(cs, "gaussian")$status, "conditional",
                 info = cs)
    expect_equal(frm_compat(cs, "poisson")$status, "conditional",
                 info = cs)
  }
  expect_match(frm_compat("exp", "gaussian")$note, "num_factor")
})

test_that("the bar-crossing refusal outranks the other grammar rules", {
  # bar_crossing x * and call_group x * have the same signature, so the
  # crossing refusal used to be lost to whichever was typed last
  expect_equal(frm_compat("bar_crossing", "call_group")$status, "refused")
  expect_equal(frm_compat("bar_crossing", "double_bar")$status, "refused")
  expect_equal(frm_compat("call_group", "double_bar")$status, "works")
})

test_that("the multivariate post-fit surface is declared as it behaves", {
  for (st in c("mvbf", "rescor")) {
    expect_equal(frm_compat(st, "fitted")$status, "refused", info = st)
    expect_equal(frm_compat(st, "predict")$status, "works", info = st)
    expect_equal(frm_compat(st, "simulate")$status, "refused", info = st)
    expect_equal(frm_compat(st, "residuals_osa")$status, "refused",
                 info = st)
    # these run: they address the outer parameter vector, not one
    # response
    expect_equal(frm_compat(st, "confint_profile")$status, "works",
                 info = st)
    expect_equal(frm_compat(st, "hypothesis_profile")$status, "works",
                 info = st)
    # frm_sample was asserted here too; the feature and its rules moved
    # to frmtmb.sample, whose own suite makes the same assertion
  }
})

# ----------------------------------------- declared refusals really refuse

test_that("declared refusals at the frame stage really refuse", {
  set.seed(2)
  n <- 60
  d <- data.frame(
    y = rnorm(n, 5), y2 = rnorm(n), x = rnorm(n),
    s = runif(n, 0.2, 0.5), w = runif(n, 0.5, 2),
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

test_that("the multivariate declarations match a multivariate fit", {
  set.seed(4)
  n <- 60
  d <- data.frame(y = rnorm(n, 2), y2 = rnorm(n, -1), x = rnorm(n))
  fit <- frm(bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), data = d,
             family = gaussian())

  expect_error(fitted(fit), "multivariate")
  expect_error(simulate(fit), "multivariate")
  expect_error(residuals(fit), "multivariate")
  expect_length(predict(fit), n)

  # the inference surface is declared to work, so it has to
  ci <- confint(fit, method = "profile", parm = "y_x")
  expect_equal(unname(ci[, "est"]),
               unname(confint(fit, method = "wald", parm = "y_x")[, "est"]),
               tolerance = 1e-6)
  expect_lt(ci[, "lwr"], ci[, "upr"])
  hp <- hypothesis(fit, "y_x = 0", method = "profile")
  expect_equal(nrow(hp), 1L)
})
