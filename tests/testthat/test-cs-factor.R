# cs() on a discrete predictor, and cs() beside the same population-level
# term.
#
# Two defects, both pre-existing through 0.64.0 and both recorded in
# dev/csfactor-findings.md with the numbers.
#
# A  A SILENT WRONG ANSWER. `ord_cs_values()` called `as.numeric()` on
#    the value, so `cs(f)` on a factor fitted the INTEGER CODES: one
#    slope over 1, 2, 3 instead of one coefficient per level. On seed
#    405 that cost 0.2872 log-likelihood units against the same model
#    written with hand-built dummy columns, and a `newdata` factor was
#    recoded from its OWN levels, so a single row `factor("c")` was read
#    as level 1 and got level a's probabilities. A character column did
#    not fit at all: `as.numeric()` on it gave NAs and the optimizer
#    died on a NaN gradient. brms 2.23.0 builds treatment-contrast
#    dummies here (`dev/csfactor-log/brms.txt`).
#
# B  `y ~ x + cs(x)` is not identified: adding a constant to the
#    population-level coefficient and subtracting it from every
#    threshold's coefficient leaves the likelihood unchanged. 0.64.0
#    fitted it and reported a standard error of 2.7e5 on all three
#    coefficients. brms builds both blocks and samples the ridge; frmtmb
#    refuses, because ML has no prior to hold it.
#
# Every tolerance here is a ratio to something the run measures: the
# log-likelihood's own magnitude, or the Monte Carlo standard error of
# the draw count actually used.

skip_on_cran()

# Seed 405 is the wt-predfix reviewer's construction, the one every
# recorded number in dev/csfactor-findings.md comes from.
csf_data <- function(seed = 405, n = 500) {
  set.seed(seed)
  x <- stats::rnorm(n)
  fc <- factor(sample(c("a", "b", "c"), n, TRUE))
  eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
  p1 <- stats::plogis(-0.3 + eff)
  p2 <- (1 - p1) * stats::plogis(0.5 - eff)
  u <- stats::runif(n)
  d <- data.frame(
    x = x, fc = fc,
    yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L))
  )
  # the same design written out by hand, brms's treatment contrasts
  d$fcb <- as.numeric(d$fc == "b")
  d$fcc <- as.numeric(d$fc == "c")
  d$fch <- as.character(d$fc)
  d$fnum <- factor(as.integer(d$fc))
  d$z <- stats::rnorm(n)
  d
}

test_that("cs() on a factor is the treatment-contrast dummies", {
  d <- csf_data()
  ff <- frm(bf(yo ~ cs(fc)) + sratio(), data = d)
  fd <- frm(bf(yo ~ cs(fcb) + cs(fcc)) + sratio(), data = d)

  # the fit is the dummy fit, not a slope on the codes. Before, the two
  # differed by 0.2872 units and the factor fit had two coefficients.
  l1 <- as.numeric(logLik(ff))
  l2 <- as.numeric(logLik(fd))
  expect_lt(abs(l1 - l2), 1e-8 * max(1, abs(l2)))
  expect_identical(nrow(fixef(ff)), nrow(fixef(fd)))
  expect_lt(max(abs(fixef(ff)[, "Estimate"] - fixef(fd)[, "Estimate"])),
            1e-6 * max(abs(fixef(fd)[, "Estimate"])))

  # brms's names: one bcs_<dummy>[k] per dummy per threshold
  expect_identical(grep("^bcs", variables(ff), value = TRUE),
                   c("bcs_fcb[1]", "bcs_fcb[2]",
                     "bcs_fcc[1]", "bcs_fcc[2]"))
  expect_identical(rownames(fixef(ff)),
                   c("Intercept[1]", "Intercept[2]", "fcb[1]", "fcb[2]",
                     "fcc[1]", "fcc[2]"))
})

test_that("a cs() factor in newdata is recoded against the FIT's levels", {
  d <- csf_data()
  ff <- frm(bf(yo ~ cs(fc)) + sratio(), data = d)
  P3 <- frm_linpred(ff, newdata = data.frame(fc = factor(c("a", "b", "c"))),
                    type = "response")

  # the defect: a one-row newdata factor carries only its own level, so
  # as.numeric() read `factor("c")` as 1 and returned level a's row
  P1 <- frm_linpred(ff, newdata = data.frame(fc = factor("c")),
                    type = "response")
  expect_lt(max(abs(P1[1L, ] - P3[3L, ])), 1e-12 * max(abs(P3)))
  # not level a's row: the bound is the run's own measured separation
  # between level a and level c, so nothing absolute is written here
  sep <- max(abs(P3[3L, ] - P3[1L, ]))
  expect_gt(max(abs(P1[1L, ] - P3[1L, ])), 0.5 * sep)

  # the same through fitted(), which is the array route
  F1 <- fitted(ff, newdata = data.frame(fc = factor("c")))
  expect_lt(max(abs(as.numeric(F1[, "Estimate", ]) - P3[3L, ])),
            1e-12 * max(abs(P3)))

  # and against the hand-written dummy model at the same row
  fd <- frm(bf(yo ~ cs(fcb) + cs(fcc)) + sratio(), data = d)
  Pd <- frm_linpred(fd, newdata = data.frame(fcb = 0, fcc = 1),
                    type = "response")
  expect_lt(max(abs(P1 - Pd)), 1e-8 * max(abs(Pd)))

  # a level the fit never saw is refused rather than recoded
  expect_error(frm_linpred(ff, newdata = data.frame(fc = factor("zz")),
                           type = "response"),
               "did not see", fixed = TRUE)
})

test_that("predict() and simulate() at a cs() factor row use its level", {
  d <- csf_data()
  ff <- frm(bf(yo ~ cs(fc)) + sratio(), data = d)
  nd <- data.frame(fc = factor("c"))
  # the reference is the HAND-DUMMY fit, not the factor fit's own
  # predictor: both halves of defect A moved together, so a factor fit
  # compared against itself agreed while both were wrong
  fd <- frm(bf(yo ~ cs(fcb) + cs(fcc)) + sratio(), data = d)
  p <- as.numeric(frm_linpred(fd, newdata = data.frame(fcb = 0, fcc = 1),
                              type = "response"))

  nsim <- 20000L
  set.seed(3)
  pr <- as.numeric(predict(ff, newdata = nd, ndraws = nsim,
                           propagate_error = FALSE))
  # the tolerance IS the draw count's own Monte Carlo standard error
  mcse <- sqrt(p * (1 - p) / nsim)
  expect_lt(max(abs(pr - p) / mcse), 5)

  set.seed(4)
  s <- simulate(ff, newdata = nd, nsim = nsim)
  sh <- as.numeric(prop.table(table(factor(unlist(s[1L, ]),
                                           levels = c("1", "2", "3")))))
  expect_lt(max(abs(sh - p) / mcse), 5)
})

test_that("cs() on a character column and on numeric-looking levels", {
  d <- csf_data()
  ref <- as.numeric(logLik(frm(bf(yo ~ cs(fc)) + sratio(), data = d)))

  # before, as.numeric() on a character column gave NAs and the fit died
  fh <- frm(bf(yo ~ cs(fch)) + sratio(), data = d)
  expect_lt(abs(as.numeric(logLik(fh)) - ref), 1e-8 * max(1, abs(ref)))
  expect_identical(grep("^bcs", variables(fh), value = TRUE),
                   c("bcs_fchb[1]", "bcs_fchb[2]",
                     "bcs_fchc[1]", "bcs_fchc[2]"))

  # a factor whose levels LOOK like the codes: before, the codes and the
  # levels agreed, so the wrong answer was invisible in the numbers
  fn <- frm(bf(yo ~ cs(fnum)) + sratio(), data = d)
  expect_lt(abs(as.numeric(logLik(fn)) - ref), 1e-8 * max(1, abs(ref)))
  expect_identical(grep("^bcs", variables(fn), value = TRUE),
                   c("bcs_fnum2[1]", "bcs_fnum2[2]",
                     "bcs_fnum3[1]", "bcs_fnum3[2]"))
})

test_that("cs() on an ORDERED factor takes its polynomial contrasts", {
  d <- csf_data()
  d$fo <- factor(as.character(d$fc), levels = c("a", "b", "c"),
                 ordered = TRUE)
  ff <- frm(bf(yo ~ cs(fo)) + sratio(), data = d)
  # model.matrix() gives an ordered factor contr.poly, and so does brms
  # (dev/csfactor-log/ordered.txt): the column names are the check that
  # the contrasts are the fit's and not a hard-coded treatment coding
  expect_identical(grep("^bcs", variables(ff), value = TRUE),
                   c("bcs_fo.L[1]", "bcs_fo.L[2]",
                     "bcs_fo.Q[1]", "bcs_fo.Q[2]"))
  # a K-1 column cs() term on a K-category response is SATURATED per
  # level, so the fitted probabilities are the level's own shares
  # whatever the contrast coding is. How close they come is the
  # OPTIMIZER's business, so the bound is the treatment-coded fit's own
  # distance from the same shares: a mishandled contrast is off by
  # tenths, not by the last digits of a converged fit.
  fd <- frm(bf(yo ~ cs(fcb) + cs(fcc)) + sratio(), data = d)
  E <- unname(prop.table(table(d$fo, d$yo), 1L))
  P <- frm_linpred(ff, newdata = data.frame(
    fo = factor(c("a", "b", "c"), levels = c("a", "b", "c"),
                ordered = TRUE)), type = "response")
  Pt <- frm_linpred(fd, newdata = data.frame(fcb = c(0, 1, 0),
                                             fcc = c(0, 0, 1)),
                    type = "response")
  expect_lt(max(abs(P - E)),
            100 * max(max(abs(Pt - E)), .Machine$double.eps))
  # the two codings are the same model space
  expect_lt(abs(as.numeric(logLik(ff)) - as.numeric(logLik(fd))),
            1e-8 * max(1, abs(as.numeric(logLik(fd)))))
})

test_that("class \"b\" priors reach a cs() factor's dummy coefficients", {
  d <- csf_data()
  ff <- frm(bf(yo ~ cs(fc)) + sratio(), data = d)

  # brms lists one class "b" row per dummy, coef = "fcb", "fcc"
  pri <- as.data.frame(default_prior(ff))
  expect_true(all(c("fcb", "fcc") %in% pri$coef[pri$class == "b"]))

  # coef = names ONE dummy: its pair shrinks, the other pair does not
  flat <- fixef(ff)[, "Estimate"]
  fp <- frm(bf(yo ~ cs(fc)) + sratio(), data = d,
            prior = set_prior("normal(0, 0.05)", class = "b",
                              coef = "fcc"))
  sh <- fixef(fp)[, "Estimate"]
  expect_lt(max(abs(sh[c("fcc[1]", "fcc[2]")])),
            0.1 * max(abs(flat[c("fcc[1]", "fcc[2]")])))
  expect_gt(max(abs(sh[c("fcb[1]", "fcb[2]")])),
            0.3 * max(abs(flat[c("fcb[1]", "fcb[2]")])))

  # class "b" with no coef covers every one of them
  fa <- frm(bf(yo ~ cs(fc)) + sratio(), data = d,
            prior = set_prior("normal(0, 0.05)", class = "b"))
  cf <- fixef(fa)[, "Estimate"]
  expect_lt(max(abs(cf[c("fcb[1]", "fcb[2]", "fcc[1]", "fcc[2]")])),
            0.15 * max(abs(flat[c("fcb[1]", "fcb[2]",
                                  "fcc[1]", "fcc[2]")])))
})

test_that("conditional_effects() on a cs() factor gives the level's row", {
  d <- csf_data()
  ff <- frm(bf(yo ~ cs(fc)) + sratio(), data = d)
  # again the hand-dummy fit is the reference, for the same reason
  fd <- frm(bf(yo ~ cs(fcb) + cs(fcc)) + sratio(), data = d)
  P3 <- frm_linpred(fd, newdata = data.frame(fcb = c(0, 1, 0),
                                             fcc = c(0, 0, 1)),
                    type = "response")
  ce <- as.data.frame(conditional_effects(ff, effects = "fc")[[1L]])
  got <- matrix(NA_real_, 3L, 3L)
  for (i in 1:3) {
    for (k in 1:3) {
      got[i, k] <- ce$estimate__[as.character(ce$fc) == c("a", "b", "c")[i] &
                                   as.character(ce$cats__) ==
                                   as.character(k)]
    }
  }
  expect_lt(max(abs(got - P3)), 1e-10 * max(abs(P3)))
})

test_that("cs() beside the same population-level term is refused", {
  d <- csf_data()
  # B, the numeric case: 0.64.0 fitted this with an se of 2.7e5
  expect_error(frm(bf(yo ~ x + cs(x)) + sratio(), data = d),
               "cs(x) is not identified", fixed = TRUE)
  expect_error(frm(bf(yo ~ x + cs(x)) + cratio(), data = d),
               "is not identified", fixed = TRUE)
  expect_error(frm(bf(yo ~ x + cs(x)) + acat(), data = d),
               "is not identified", fixed = TRUE)
  # and the factor case, where the redundancy is one dummy at a time
  expect_error(frm(bf(yo ~ fc + cs(fc)) + sratio(), data = d),
               "cs(fc) is not identified", fixed = TRUE)
  expect_error(frm(bf(yo ~ fc + cs(fc)) + sratio(), data = d),
               "its column 'fcb'", fixed = TRUE)
  # a constant cs() column has nothing to separate it from a threshold
  d$one <- 1
  expect_error(frm(bf(yo ~ cs(one)) + sratio(), data = d),
               "is constant over the rows", fixed = TRUE)
  # the message names the edit that fixes it
  expect_error(frm(bf(yo ~ x + cs(x)) + sratio(), data = d),
               "Write 'x + cs(x)' as 'cs(x)' alone", fixed = TRUE)

  # 'cs(x)' alone fits the same set of distributions the refused model
  # spans, which is why the message says to write it that way
  f1 <- frm(bf(yo ~ cs(x)) + sratio(), data = d)
  expect_true(is.finite(as.numeric(logLik(f1))))
})

test_that("the identifiability check does not fire on a separable model", {
  d <- csf_data()
  # the shapes the suite and the vignettes use: a different variable on
  # each side, a factor beside a numeric cs(), two cs() terms, and a cs()
  # whose variable is in no other term
  expect_no_error(frm(bf(yo ~ z + cs(x)) + sratio(), data = d))
  expect_no_error(frm(bf(yo ~ fc + cs(x)) + sratio(), data = d))
  expect_no_error(frm(bf(yo ~ x + cs(fc)) + sratio(), data = d))
  expect_no_error(frm(bf(yo ~ cs(fc) + cs(x)) + sratio(), data = d))
  expect_no_error(frm(bf(yo ~ cs(x)) + sratio(), data = d))
  # an interaction of two cs() columns is NOT in their span
  expect_no_error(frm(bf(yo ~ cs(x) + cs(x:z)) + sratio(), data = d))
})

test_that("cs() works in a multivariate fit, in every response", {
  # `?frm` said cs() was "not available ... in a multivariate fit". It
  # was never true, and this pins the page to the behavior.
  d <- csf_data()
  set.seed(778)
  d$yo2 <- sample(1:3, nrow(d), TRUE)
  mv <- frm(bf(yo ~ x + cs(fc)) + bf(yo2 ~ z + cs(fc)) + sratio(), data = d)
  expect_identical(grep("^bcs", variables(mv), value = TRUE),
                   c("bcs_yo_fcb[1]", "bcs_yo_fcb[2]",
                     "bcs_yo_fcc[1]", "bcs_yo_fcc[2]",
                     "bcs_yo2_fcb[1]", "bcs_yo2_fcb[2]",
                     "bcs_yo2_fcc[1]", "bcs_yo2_fcc[2]"))
  expect_true(all(c("yo_fcb[1]", "yo_fcc[2]", "yo2_fcb[1]",
                    "yo2_fcc[2]") %in% rownames(fixef(mv))))

  # with no rescor the multivariate likelihood IS the product of the two,
  # so this is an IDENTITY and the residual is the optimizer's noise
  u1 <- frm(bf(yo ~ x + cs(fc)) + sratio(), data = d)
  u2 <- frm(bf(yo2 ~ z + cs(fc)) + sratio(), data = d)
  tot <- as.numeric(logLik(u1)) + as.numeric(logLik(u2))
  expect_lt(abs(as.numeric(logLik(mv)) - tot), 1e-6 * max(1, abs(tot)))

  # newdata prediction reaches one response's cs() columns
  nd <- data.frame(x = 0, z = 0, fc = factor("c", levels = c("a", "b", "c")))
  P <- frm_linpred(mv, newdata = nd, resp = "yo", type = "response")
  expect_equal(unname(rowSums(P)), 1, tolerance = 1e-12)

  # and ONE bad predictor is still refused, naming which one
  expect_error(frm(bf(yo ~ x + cs(x)) + bf(yo2 ~ z + cs(fc)) + sratio(),
                   data = d), "'yo.mu'", fixed = TRUE)
})

test_that("mo() beside cs() on the same predictor is refused", {
  # The rank test cannot see the mo() column, which is a zero placeholder
  # at assembly, so mo() gets a test against the code-indicator basis
  # instead. Unfixed, `yo ~ mo(m) + cs(m)` fitted with 3 extra degrees of
  # freedom buying 6.4e-09 of log likelihood and 9 of 9 NaN standard
  # errors (dev/csfactor-rev-log/rev-gap-lane.log, seed 1907).
  set.seed(1907)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$m <- factor(sample(1:4, n, TRUE), ordered = TRUE)
  d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
  d$yo <- sample(1:3, n, TRUE)
  expect_error(frm(bf(yo ~ mo(m) + cs(m)) + sratio(), data = d),
               "not identified together", fixed = TRUE)
  expect_error(frm(bf(yo ~ mo(m) + cs(m)) + sratio(), data = d),
               "whatever its simplex", fixed = TRUE)

  # and it must NOT fire where the mo() direction is still free. These
  # are the constructions the check would fire on if it used the whole
  # indicator basis of `m` as the mo() column rather than the span every
  # simplex reaches.
  expect_no_error(frm(bf(yo ~ mo(m)) + sratio(), data = d))
  expect_no_error(frm(bf(yo ~ cs(m)) + sratio(), data = d))
  expect_no_error(frm(bf(yo ~ mo(m) + cs(x)) + sratio(), data = d))
  expect_no_error(frm(bf(yo ~ mo(m) + cs(f)) + sratio(), data = d))
  # a COARSENING of m: cs(mc) spans only 2 of m's 4 values, so the
  # monotonic shape is still identified
  d$mc <- factor(ifelse(as.integer(d$m) <= 2L, "lo", "hi"))
  expect_no_error(frm(bf(yo ~ mo(m) + cs(mc)) + sratio(), data = d))
  # mo(m):z is z times a function of m; cs(m) cannot absorb it, so the
  # INTERACTION survives and the message says which spelling to use.
  # `mo(m) * z` expands to the main effect plus the interaction, and it is
  # the main effect that is refused.
  expect_no_error(frm(bf(yo ~ mo(m):z + cs(m)) + sratio(), data = d))
  expect_no_error(frm(bf(yo ~ z + mo(m):z + cs(m)) + sratio(), data = d))
  expect_error(frm(bf(yo ~ mo(m) * z + cs(m)) + sratio(), data = d),
               "not identified together", fixed = TRUE)
  expect_error(frm(bf(yo ~ mo(m) * z + cs(m)) + sratio(), data = d),
               "write 'z + mo(m):z + cs(m)'", fixed = TRUE)
})

test_that("the refusal names the columns and the edit that works", {
  set.seed(1907)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$yo <- sample(1:3, n, TRUE)
  # a plain column: delete the term
  expect_error(frm(bf(yo ~ x + cs(x)) + sratio(), data = d),
               "predictor ('x')", fixed = TRUE)
  expect_error(frm(bf(yo ~ x + cs(x)) + sratio(), data = d),
               "Write 'x + cs(x)' as 'cs(x)' alone", fixed = TRUE)
  # a polynomial basis: there is no `x` term to delete, and the advice
  # must not tell the reader to delete one
  e <- tryCatch(frm(bf(yo ~ poly(x, 2) + cs(x)) + sratio(), data = d),
                error = conditionMessage)
  expect_match(e, "poly(x, 2)1", fixed = TRUE)
  expect_match(e, "I(x^2) + cs(x)", fixed = TRUE)
  expect_false(grepl("Write 'x + cs(x)'", e, fixed = TRUE))
  # a smooth: no spelling keeps both
  es <- tryCatch(frm(bf(yo ~ s(x) + cs(x)) + sratio(), data = d),
                 error = conditionMessage)
  expect_match(es, "SMOOTH's unpenalized linear part", fixed = TRUE)
  expect_false(grepl("I(x^2)", es, fixed = TRUE))
  expect_match(e, "POLYNOMIAL basis", fixed = TRUE)
  # two population-level columns spanning one cs() column: name both, and
  # give GENERIC advice. The poly worked example used to be pasted here
  # too, advertising a poly() this formula does not contain.
  e2 <- tryCatch(frm(bf(yo ~ x + z + cs(I(x + z))) + sratio(), data = d),
                 error = conditionMessage)
  expect_match(e2, "('x', 'z')", fixed = TRUE)
  expect_match(e2, "nothing to delete by that name", fixed = TRUE)
  expect_match(e2, "drop either the cs() term or the population-level",
               fixed = TRUE)
  expect_false(grepl("poly", e2, fixed = TRUE))
  expect_false(grepl("I(x^2)", e2, fixed = TRUE))
  # each worked example stays in its own branch
  expect_false(grepl("nothing to delete by that name", e, fixed = TRUE))
  expect_false(grepl("POLYNOMIAL", es, fixed = TRUE))
})

test_that("an unused cs() factor level drops its column, as in brms", {
  # assemble_frame() builds its model frame with drop.unused.levels =
  # TRUE, upstream of the cs() design, so an absent level never becomes
  # an all-zero dummy and the constant-column refusal never sees one.
  # brms drops it too (Kcs = 1, Xcs = fub).
  set.seed(405)
  n <- 240
  d <- data.frame(x = stats::rnorm(n))
  d$f <- factor(c(rep("a", 114), rep("b", 126)),
                levels = c("a", "b", "c"))
  d$yo <- sample(1:3, n, TRUE)
  expect_identical(nlevels(d$f), 3L)          # the level IS declared
  ff <- frm(bf(yo ~ x + cs(f)) + sratio(), data = d)
  expect_length(ff$frame$linpreds[["yo.mu"]][["cs"]], 1L)
  expect_identical(grep("^bcs", variables(ff), value = TRUE),
                   c("bcs_fb[1]", "bcs_fb[2]"))
})

test_that("influence() warns when a deletion drops a cs() coefficient", {
  # The refit SUCCEEDS with one column fewer, so the failed-refit count
  # cannot see it: those cells stay NA and cooks.distance() is NA for
  # that unit. It used to do that silently.
  set.seed(90291)
  n <- 120
  d <- data.frame(x = stats::rnorm(n))
  d$f <- factor(c(rep("a", 60), rep("b", 59), "c"),
                levels = c("a", "b", "c"))
  d$yo <- sample(1:3, n, TRUE)
  ff <- frm(bf(yo ~ x + cs(f)) + sratio(), data = d)
  expect_warning(inf <- influence(ff, force = TRUE),
                 "cooks.distance() is NA", fixed = TRUE)
  cd <- suppressWarnings(cooks.distance(influence(ff, force = TRUE)))
  # exactly the unit whose deletion removes the only row of level c
  expect_identical(unname(which(is.na(cd))), 120L)
  # the warning names the level, not only the internal parameter
  w <- tryCatch(influence(ff, force = TRUE),
                warning = conditionMessage)
  expect_match(w, "cs fc", fixed = TRUE)
  expect_match(w, "'120'", fixed = TRUE)
})

test_that("a cs() factor survives a re-fit through the stored frame", {
  # the model matrix and the levels travel with the frame, so a path
  # that rebuilds the objective from it sees the same columns
  d <- csf_data()
  ff <- frm(bf(yo ~ cs(fc)) + sratio(), data = d)
  fr <- frm(bf(yo ~ cs(fc)) + sratio(), data = d, dry_run = "frame")
  cs <- fr$linpreds[["yo.mu"]][["cs"]]
  expect_length(cs, 2L)
  expect_identical(vapply(cs, `[[`, "", "label"), c("csfcb", "csfcc"))
  # a 0/1 dummy is exact in both places, so this is bitwise on purpose
  expect_identical(unname(cs[[1L]][["vals"]]), as.numeric(d$fc == "b"))
  expect_identical(unname(cs[[2L]][["vals"]]), as.numeric(d$fc == "c"))
  mm <- fr$linpreds[["yo.mu"]][["cs_mm"]]
  expect_length(mm, 1L)
  expect_identical(mm[[1L]][["xlevels"]][["fc"]], c("a", "b", "c"))
  expect_identical(mm[[1L]][["colnames"]], c("fcb", "fcc"))
  nll <- frmtmb:::build_objective(ff$frame)
  expect_lt(abs(-nll(ff$estimates) - as.numeric(logLik(ff))),
            1e-8 * max(1, abs(as.numeric(logLik(ff)))))
})
