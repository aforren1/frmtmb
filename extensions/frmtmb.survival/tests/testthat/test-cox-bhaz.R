# brms's bhaz() addition term. brms 2.23.0 writes the basis of the
# cox() baseline on the response side, time | bhaz(df = 6) ~ x, and
# this package takes it as an option of cox() through
# frmtmb_register_aterm(option =). Each spelling must give the fit of
# the cox() argument of the same name, and what brms takes and this
# family cannot fit, bhaz(gr = ), must be refused by name.

bhaz_data <- function(seed = 11, n = 200) {
  set.seed(seed)
  x <- stats::rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
             g = factor(rep(c("a", "b", "c", "d"), length.out = n)),
             z = stats::rnorm(n))
}

# the two spellings are one model: the same objective, so the same fit
# to the last bit
expect_same_fit <- function(a, b) {
  expect_identical(logLik(a), logLik(b))
  expect_identical(a$estimates, b$estimates)
  expect_identical(vcov(a), vcov(b))
}

test_that("bhaz() gives the fit of the cox() argument of the same name", {
  dd <- bhaz_data()
  pairs <- list(
    list(time | cens(cens) + bhaz(df = 6) ~ x, cox(df = 6)),
    list(time | cens(cens) + bhaz(knots = c(0.5, 1.1),
                                  Boundary.knots = c(0, 16)) ~ x,
         cox(knots = c(0.5, 1.1), Boundary.knots = c(0, 16))),
    list(time | cens(cens) + bhaz(Boundary.knots = c(0, 16)) ~ x,
         cox(Boundary.knots = c(0, 16))),
    list(time | cens(cens) + bhaz(df = 5, degree = 2,
                                  intercept = FALSE) ~ x,
         cox(df = 5, degree = 2, intercept = FALSE)))
  for (p in pairs) {
    a <- frm(bf(p[[1L]]), family = cox(), data = dd)
    b <- frm(bf(time | cens(cens) ~ x), family = p[[2L]], data = dd)
    expect_same_fit(a, b)
  }
  # and the basis is what was asked for, not the default one
  a <- frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(),
           data = dd)
  expect_length(cox_baseline(a), 6L)
  k <- frm(bf(time | cens(cens) + bhaz(knots = c(0.5, 1.1)) ~ x),
           family = cox(), data = dd, dry_run = "frame")
  sp <- k$aterm_values$time$bhaz_spec
  expect_identical(sp$internal, c(0.5, 1.1))
  expect_identical(ncol(k$aterm_values$time$Zbhaz), 6L)
})

test_that("without bhaz(), or with only its defaults, cox() is unchanged", {
  # the absent case: the term must change nothing it does not name
  dd <- bhaz_data(seed = 12)
  ref <- frm(bf(time | cens(cens) ~ x), family = cox(), data = dd)
  expect_same_fit(frm(bf(time | cens(cens) + bhaz() ~ x), family = cox(),
                      data = dd), ref)
  expect_same_fit(frm(bf(time | cens(cens) + bhaz(gr = NA) ~ x),
                      family = cox(), data = dd), ref)
  expect_same_fit(frm(bf(time | cens(cens) + bhaz(df = 5) ~ x),
                      family = cox(), data = dd), ref)
  expect_length(cox_baseline(ref), 5L)
  # a cox() argument bhaz() does not name is kept
  r6 <- frm(bf(time | cens(cens) ~ x), family = cox(df = 6), data = dd)
  expect_same_fit(frm(bf(time | cens(cens) + bhaz() ~ x),
                      family = cox(df = 6), data = dd), r6)
  # and one it names with the same value agrees
  expect_same_fit(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x),
                      family = cox(df = 6), data = dd), r6)
  # the string and brmsfamily() spellings of the family take it too
  expect_same_fit(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x),
                      family = "cox", data = dd), r6)
  expect_same_fit(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x),
                      family = brmsfamily("cox"), data = dd), r6)
})

test_that("bhaz() reads its arguments in the formula's environment", {
  dd <- bhaz_data(seed = 13)
  kn <- c(0.5, 1.1)
  a <- frm(bf(time | cens(cens) + bhaz(knots = kn) ~ x), family = cox(),
           data = dd)
  b <- frm(bf(time | cens(cens) ~ x), family = cox(knots = c(0.5, 1.1)),
           data = dd)
  expect_same_fit(a, b)
  # the knots are not a data column
  expect_false("kn" %in% a$frame$data_vars)
  expect_error(frm(bf(time | cens(cens) + bhaz(df = nope) ~ x),
                   family = cox(), data = dd),
               "bhaz(): cannot evaluate df = nope", fixed = TRUE)
})

test_that("a bhaz() fit survives update() and new data", {
  dd <- bhaz_data(seed = 14)
  a <- frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(),
           data = dd)
  up <- update(a, . ~ . + z)
  ref <- frm(bf(time | cens(cens) ~ x + z), family = cox(df = 6),
             data = dd)
  expect_same_fit(up, ref)
  nd <- dd[1:5, ]
  expect_identical(frm_linpred(a, newdata = nd, type = "link"),
                   frm_linpred(frm(bf(time | cens(cens) ~ x),
                                   family = cox(df = 6), data = dd),
                               newdata = nd, type = "link"))
})

test_that("bhaz() and cox() that disagree are refused by name", {
  dd <- bhaz_data(seed = 15)
  fo <- function(lhs) stats::as.formula(paste(lhs, "~ x"))
  expect_error(frm(bf(fo("time | cens(cens) + bhaz(df = 6)")),
                   family = cox(df = 5), data = dd),
               "bhaz(df = 6) and cox(df = 5) set the baseline basis twice",
               fixed = TRUE)
  expect_error(frm(bf(fo("time | bhaz(degree = 2)")),
                   family = cox(degree = 3), data = dd),
               "bhaz(degree = 2) and cox(degree = 3)", fixed = TRUE)
  expect_error(frm(bf(fo("time | bhaz(knots = c(0.5, 1.1))")),
                   family = cox(knots = c(0.5, 1.2)), data = dd),
               "bhaz(knots = c(0.5, 1.1)) and cox(knots = c(0.5, 1.2))",
               fixed = TRUE)
  # df and knots size one basis, so the two must agree wherever each
  # is written; brms ignores df when knots is given
  expect_error(frm(bf(fo("time | bhaz(df = 8, knots = c(0.5, 1.1))")),
                   family = cox(), data = dd),
               "bhaz(df =) is 8, and bhaz(knots =) gives 2 internal knots",
               fixed = TRUE)
  expect_error(frm(bf(fo("time | bhaz(knots = c(0.5, 1.1))")),
                   family = cox(df = 4), data = dd),
               "cox(df =) is 4, and bhaz(knots =) gives 2", fixed = TRUE)
  expect_error(cox(df = 4, knots = c(0.5, 1.1)),
               "cox(df =) is 4, and cox(knots =) gives 2", fixed = TRUE)
  # the absent case: knots that agree with df are taken
  expect_same_fit(frm(bf(fo("time | bhaz(df = 6, knots = c(0.5, 1.1))")),
                      family = cox(), data = dd),
                  frm(bf(time ~ x), family = cox(knots = c(0.5, 1.1)),
                      data = dd))
})

test_that("bhaz(gr = ), a stratified baseline, is refused by name", {
  dd <- bhaz_data(seed = 16)
  for (lhs in c("time | cens(cens) + bhaz(gr = g)", "time | bhaz(g)",
                "time | bhaz(gr = g, df = 4)")) {
    expect_error(frm(bf(stats::as.formula(paste(lhs, "~ x"))),
                     family = cox(), data = dd),
                 "a baseline hazard of its own (a stratified Cox model)",
                 fixed = TRUE)
  }
  e <- tryCatch(frm(bf(time | bhaz(gr = g) ~ x), family = cox(),
                    data = dd), error = function(e) e)
  # it says what to do instead
  expect_match(conditionMessage(e), "Fit each stratum as a model of its own",
               fixed = TRUE)
  expect_match(conditionMessage(e), "or add g to the formula", fixed = TRUE)
  expect_s3_class(e, "frmtmb_survival_error")
})

test_that("bhaz() refuses what it cannot fit, by name", {
  dd <- bhaz_data(seed = 17)
  fo <- function(lhs) bf(stats::as.formula(paste(lhs, "~ x")))
  expect_error(frm(fo("time | bhaz(foo = 1)"), family = cox(), data = dd),
               "bhaz() takes brms's arguments for the baseline basis",
               fixed = TRUE)
  # brms extrapolates past boundary knots inside the data, to a hazard
  # that can be negative
  expect_error(frm(fo("time | bhaz(Boundary.knots = c(0.5, 10))"),
                   family = cox(), data = dd),
               "must cover every observed time", fixed = TRUE)
  expect_error(frm(fo("time | bhaz(knots = c(0.5, 20))"), family = cox(),
                   data = dd),
               "must lie strictly inside its boundary knots", fixed = TRUE)
  expect_error(cox(knots = c(0.5, 20), Boundary.knots = c(0, 16)),
               "must lie strictly inside its boundary knots", fixed = TRUE)
  expect_error(cox(knots = c(1, 1)), "takes distinct finite numbers",
               fixed = TRUE)
  expect_error(cox(Boundary.knots = c(2, 1)),
               "takes two finite numbers", fixed = TRUE)
  expect_error(frm(fo("time | bhaz(df = 3)"), family = cox(), data = dd),
               "bhaz(): df must be at least degree + 1 (4)", fixed = TRUE)
  # cox() alone is refused before the fit, at frame assembly; the
  # constructor leaves it to the merge with a bhaz() that may follow
  expect_error(frm(fo("time"), family = cox(df = 3), data = dd,
                   dry_run = "frame"),
               "cox(): df must be at least degree + 1 (4)", fixed = TRUE)
  expect_error(frm(fo("time"), family = cox(df = 2, intercept = FALSE),
                   data = dd, dry_run = "frame"),
               "cox(): df must be at least degree (3) without an intercept",
               fixed = TRUE)
  expect_error(frm(fo("time | bhaz(degree = 4)"), family = cox(df = 4),
                   data = dd, dry_run = "frame"),
               "bhaz(): df must be at least degree + 1 (5)", fixed = TRUE)
  # brms truncates a fractional df; frmtmb refuses it
  expect_error(frm(fo("time | bhaz(df = 4.5)"), family = cox(), data = dd),
               "`bhaz(df =)` must be a single whole number", fixed = TRUE)
  expect_error(frm(fo("time | bhaz(df = 6) + bhaz(df = 6)"),
                   family = cox(), data = dd),
               "Duplicated addition term `bhaz()`", fixed = TRUE)
  # brms takes bhaz() for cox() alone
  expect_error(frm(fo("time | bhaz(df = 6)"), family = gaussian(),
                   data = dd),
               "Addition term `bhaz()` is not supported for family 'gaussian'",
               fixed = TRUE)
  # the absent case of the boundary check: knots that cover the data,
  # and knots exactly at its ends
  bk <- c(0, 2 * max(dd$time))
  expect_no_error(frm(fo("time | bhaz(Boundary.knots = bk)"),
                      family = cox(), data = dd, dry_run = "frame"))
  bk <- range(dd$time)
  expect_no_error(frm(fo("time | bhaz(Boundary.knots = bk)"),
                      family = cox(), data = dd, dry_run = "frame"))
})

test_that("bhaz() can make a cox() df valid: the basis is checked merged", {
  dd <- bhaz_data(seed = 19)
  # df 3 is too few with an intercept and enough without one
  a <- frm(bf(time | cens(cens) + bhaz(intercept = FALSE) ~ x),
           family = cox(df = 3), data = dd)
  b <- frm(bf(time | cens(cens) ~ x),
           family = cox(df = 3, intercept = FALSE), data = dd)
  expect_same_fit(a, b)
  expect_length(cox_baseline(a), 3L)
  a2 <- frm(bf(time | cens(cens) + bhaz(degree = 2) ~ x),
            family = cox(df = 3), data = dd)
  expect_same_fit(a2, frm(bf(time | cens(cens) ~ x),
                          family = cox(df = 3, degree = 2), data = dd))
})

test_that("bhaz() reads argument names as brms 2.23.0 reads them", {
  # dev/bhaz-brms-partial.R: what brms takes, frmtmb takes; what brms
  # drops without a word, frmtmb refuses and names the argument meant
  dd <- bhaz_data(seed = 20)
  fo <- function(lhs) bf(stats::as.formula(paste(lhs, "~ x")))
  fit <- function(lhs, fam = cox()) frm(fo(lhs), family = fam, data = dd)
  ref6 <- fit("time | cens(cens)", cox(df = 6))
  # brms: d = 6 is df = 6 (resp_bhaz() matches it), kn and deg reach
  # splines2's knots and degree, g is gr
  expect_same_fit(fit("time | cens(cens) + bhaz(d = 6)"), ref6)
  expect_same_fit(fit("time | cens(cens) + bhaz(kn = c(0.5, 1.1))"),
                  fit("time | cens(cens)", cox(knots = c(0.5, 1.1))))
  expect_same_fit(fit("time | cens(cens) + bhaz(k = c(0.5, 1.1))"),
                  fit("time | cens(cens)", cox(knots = c(0.5, 1.1))))
  expect_same_fit(fit("time | cens(cens) + bhaz(deg = 2)"),
                  fit("time | cens(cens)", cox(degree = 2)))
  expect_same_fit(fit("time | cens(cens) + bhaz(NA, 6)"), ref6)
  expect_error(fit("time | bhaz(g = g)"), "a stratified Cox model",
               fixed = TRUE)
  # warn.outside reaches splines2 in brms and changes no basis
  ref <- fit("time | cens(cens)")
  expect_same_fit(fit("time | cens(cens) + bhaz(warn.outside = FALSE)"),
                  ref)
  expect_same_fit(fit("time | cens(cens) + bhaz(warn = FALSE)"), ref)
  # brms keeps intercept = TRUE and the default boundary knots for these,
  # a different model from the one written
  for (nm in c("int", "inter", "i")) {
    expect_error(fit(paste0("time | bhaz(", nm, " = FALSE)")),
                 paste0("bhaz(", nm, " = ) is not an argument of bhaz(): ",
                        "write intercept ="), fixed = TRUE)
  }
  for (nm in c("Boundary", "Bound", "B")) {
    expect_error(fit(paste0("time | bhaz(", nm, " = c(0, 30))")),
                 paste0("bhaz(", nm, " = ) is not an argument of bhaz(): ",
                        "write Boundary.knots ="), fixed = TRUE)
  }
  # splines2 arguments that change brms's basis, and cox() cannot take
  for (a in c("periodic = TRUE", "per = TRUE", "derivs = 1",
              "integral = TRUE")) {
    expect_error(fit(paste0("time | bhaz(", a, ")")),
                 "changes the baseline basis in brms", fixed = TRUE)
  }
  # ambiguous in brms too, and names brms drops
  expect_error(fit("time | bhaz(de = 2)"), "could be any of degree or derivs",
               fixed = TRUE)
  expect_error(fit("time | bhaz(Knots = c(0.5, 1.1))"),
               "`Knots =` is none of them", fixed = TRUE)
  expect_error(fit("time | bhaz(df = 6, d = 7)"),
               "could be any of degree or derivs", fixed = TRUE)
  expect_error(fit("time | bhaz(df = 6, df = 7)"),
               paste("bhaz() sets df twice, as `df =` and earlier in the",
                     "same call. brms stops on this too"), fixed = TRUE)
  # a full name and its shortened form: brms keeps the full one
  # (dev/bhaz-brms-partial.R: degree 2, knots 0.5 1.1); frmtmb refuses
  expect_error(fit("time | bhaz(degree = 2, deg = 3)"),
               "brms keeps `degree =` and drops `deg =` without a message",
               fixed = TRUE)
  expect_error(fit("time | bhaz(knots = c(0.5, 1.1), kn = c(1, 2))"),
               "brms keeps `knots =` and drops `kn =` without a message",
               fixed = TRUE)
  # brms ignores a third unnamed argument
  expect_error(fit("time | bhaz(NA, 6, c(0.5, 1.1))"),
               "at most two unnamed arguments", fixed = TRUE)
  # bhaz(6) is gr = 6 in brms, one stratum named 6
  e <- tryCatch(fit("time | bhaz(6)"), error = function(e) e)
  expect_match(conditionMessage(e), "the first argument of bhaz() is gr",
               fixed = TRUE)
  expect_match(conditionMessage(e), "Write bhaz(df = 6)", fixed = TRUE)
  expect_false(grepl("add 6 to the formula", conditionMessage(e),
                     fixed = TRUE))
})

test_that("frm_compat() says what frm() does with bhaz()", {
  rows <- frm_compat("bhaz()")
  st <- function(b) rows$status[rows$feature_b == b]
  expect_identical(st("cox"), "works")
  expect_identical(st("royston_parmar"), "refused")
  expect_identical(st("gaussian"), "refused")
  dd <- bhaz_data(seed = 21)
  expect_error(frm(bf(time | bhaz(df = 6) ~ x), family = royston_parmar(),
                   data = dd),
               "family 'royston_parmar' has none", fixed = TRUE)
})

test_that("cox(bhaz = ), brms's spelling before 2.22.0, is refused", {
  expect_error(cox(bhaz = list(df = 6)),
               "cox(bhaz = ) is the spelling of brms before 2.22.0",
               fixed = TRUE)
  expect_error(cox(bhaz = list(df = 6)), "time | bhaz(df = 6) ~ x",
               fixed = TRUE)
})

test_that("the bhaz() basis equals brms 2.23.0's standata()", {
  skip_if_not_installed("brms", "2.23.0")
  skip_if_not_installed("splines2")
  dd <- bhaz_data(seed = 18)
  lhs <- c("time | cens(cens) + bhaz(df = 6)",
           paste("time | cens(cens) + bhaz(knots = c(0.5, 1.1),",
                 "Boundary.knots = c(0, 16))"),
           "time | cens(cens) + bhaz(df = 5, degree = 2, intercept = FALSE)",
           "time | cens(cens)")
  for (l in lhs) {
    fo <- stats::as.formula(paste(l, "~ x"))
    sd <- brms::standata(brms::bf(fo), data = dd, family = brms::cox())
    fr <- frm(bf(fo), family = cox(), data = dd, dry_run = "frame")
    av <- fr$aterm_values$time
    for (nm in c("Zbhaz", "Zcbhaz")) {
      ref <- matrix(as.numeric(sd[[nm]]), nrow(sd[[nm]]))
      expect_identical(dim(av[[nm]]), dim(ref), info = l)
      # a few ulps of the basis's own scale: two implementations of one
      # recursion
      expect_lte(max(abs(av[[nm]] - ref)),
                 16 * .Machine$double.eps * max(abs(ref)))
    }
  }
})
