# Nested and multi-membership group labels, as brms 2.23.0 writes them:
# (1 + x | g1/g2) has the groups g1 and g1:g2, in that order, and
# (1 | mm(m1, m2)) the group mmm1m2. The expected names are brms's own
# from a fit of the same model (dev/prilabel-brms-fit.R,
# dev/prilabel-log/brms-fit.txt and brms-vars.txt).

lab_data <- function(n = 300, seed = 5) {
  set.seed(seed)
  g1 <- factor(sample(sprintf("p%d", 1:6), n, TRUE))
  g2 <- factor(sample(sprintf("q%d", 1:5), n, TRUE))
  m1 <- factor(sample(sprintf("m%02d", 1:10), n, TRUE))
  m2 <- factor(sample(sprintf("m%02d", 1:10), n, TRUE))
  g <- factor(sample(sprintf("g%02d", 1:12), n, TRUE))
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(6, 0, 0.7)[as.integer(g1)] +
    rnorm(30, 0, 0.5)[(as.integer(g1) - 1L) * 5L + as.integer(g2)] +
    0.5 * (rnorm(10, 0, 0.5)[as.integer(m1)] +
             rnorm(10, 0, 0.5)[as.integer(m2)]) + rnorm(n, 0, 0.6)
  data.frame(y, x, g1, g2, m1, m2, g)
}

lab_fit <- function(dd, ...) {
  frm(bf(y ~ x + (1 + x | g1/g2) + (1 | mm(m1, m2)), sigma ~ (1 | g)),
      family = gaussian(), data = dd, ...)
}

test_that("nested and mm groups carry brms's labels everywhere", {
  dd <- lab_data()
  fit <- lab_fit(dd)
  # brms's variables() of the group-level parameters: every sd, then
  # every cor
  v <- variables(fit)
  expect_identical(
    v[grepl("^(sd|cor)_", v)],
    c("sd_g1__Intercept", "sd_g1__x", "sd_g1:g2__Intercept", "sd_g1:g2__x",
      "sd_mmm1m2__Intercept", "sd_g__sigma_Intercept",
      "cor_g1__Intercept__x", "cor_g1:g2__Intercept__x"))
  grp <- c("g1", "g1:g2", "mmm1m2", "g")
  expect_identical(names(ranef(fit)), grp)
  expect_identical(names(VarCorr(fit)), grp)
  expect_identical(names(ngrps(fit)), grp)
  expect_identical(names(summary(fit)$random), grp)
  # brms's levels of an interaction group join its parts with "_"
  expect_identical(head(rownames(ranef(fit)[["g1:g2"]]), 2),
                   c("p1_q1", "p1_q2"))
  expect_identical(rownames(coef(fit)$mu[["g1:g2"]]),
                   rownames(ranef(fit)[["g1:g2"]]))
  p <- as.data.frame(default_prior(
    bf(y ~ x + (1 + x | g1/g2) + (1 | mm(m1, m2))), data = dd,
    family = gaussian()))
  expect_identical(unique(p$group[p$class == "sd" & nzchar(p$group)]),
                   c("g1", "g1:g2", "mmm1m2"))
  # a new level is named by brms's group, and spelled as brms spells it
  nd <- dd[1:2, ]
  nd$g2 <- factor("qNEW")
  expect_error(fitted(fit, newdata = nd),
               paste0("grouping factor `g1:g2`: ", dd$g1[1], "_qNEW"),
               fixed = TRUE)
  nd <- dd[1:2, ]
  nd$m1 <- factor("mNEW")
  expect_error(fitted(fit, newdata = nd),
               "multi-membership factor `mmm1m2`", fixed = TRUE)
})

test_that("set_prior() takes brms's group labels and only those", {
  dd <- lab_data()
  pr <- set_prior("normal(0, 0.5)", class = "sd", group = "g1:g2") +
    set_prior("normal(0, 0.3)", class = "sd", group = "mmm1m2")
  fit <- lab_fit(dd, prior = pr)
  tab <- as.data.frame(prior_summary(fit))
  row <- function(g) tab$prior[tab$class == "sd" & tab$group == g &
                                 !nzchar(tab$coef)]
  expect_identical(row("g1:g2"), "normal(0, 0.5)")
  expect_identical(row("mmm1m2"), "normal(0, 0.3)")
  # brms refuses reformulas's order and the term as written
  # (dev/prilabel-log/brms-ref2.txt); no alias is kept
  for (g in c("g2:g1", "mm(m1, m2)")) {
    expect_error(lab_fit(dd, prior = set_prior("normal(0, 1)", class = "sd",
                                               group = g)),
                 "No random-effect SDs match", fixed = TRUE)
  }
  # an interaction written g2:g1 is brms's g2:g1, as written
  p <- as.data.frame(default_prior(bf(y ~ x + (1 | g2:g1)), data = dd,
                                   family = gaussian()))
  expect_identical(unique(p$group[p$class == "sd" & nzchar(p$group)]),
                   "g2:g1")
})

test_that("g1/g2 is the model (1 | g1) + (1 | g1:g2)", {
  dd <- lab_data()
  a <- frm(bf(y ~ x + (1 | g1/g2)), family = gaussian(), data = dd)
  b <- frm(bf(y ~ x + (1 | g1) + (1 | g1:g2)), family = gaussian(),
           data = dd)
  expect_identical(variables(a), variables(b))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)))
  # three levels deep, as brms: g1, g1:g2, g1:g2:g3
  dd$g3 <- factor(sample(c("r1", "r2"), nrow(dd), TRUE))
  c3 <- frm(bf(y ~ x + (1 | g1/g2/g3)), family = gaussian(), data = dd)
  expect_identical(names(VarCorr(c3)),
                   c("g1", "g1:g2", "g1:g2:g3", "residual__"))
  # a covariance wrapper and || keep their structure on every piece
  d2 <- frm(bf(y ~ x + (1 + x || g1/g2)), family = gaussian(), data = dd)
  expect_identical(grep("^cor_", variables(d2), value = TRUE),
                   character(0))
  expect_identical(names(VarCorr(d2)), c("g1", "g1:g2", "residual__"))
  # brms refuses g1/g2 written twice over one coefficient
  expect_error(frm(bf(y ~ x + (1 | g1/g2) + (1 | g1:g2)),
                   family = gaussian(), data = dd),
               "Duplicated group-level effects", fixed = TRUE)
})

test_that("gr(g, id = ) is brms's |ID| key", {
  dd <- lab_data()
  a <- frm(bf(y ~ x + (1 | gr(g1, id = "a")) + (0 + x | gr(g1, id = "a"))),
           family = gaussian(), data = dd)
  b <- frm(bf(y ~ x + (1 | a | g1) + (0 + x | a | g1)), family = gaussian(),
           data = dd)
  # one correlated block over both terms, as brms builds it
  # (dev/prilabel-log/brms-ref.txt: one cor row, group g1)
  expect_identical(variables(a), variables(b))
  expect_true("cor_g1__Intercept__x" %in% variables(a))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)))
  # brms's id = NA is no id: two blocks, no correlation
  c0 <- frm(bf(y ~ x + (1 | gr(g1, id = NA)) + (0 + x | gr(g1, id = NA))),
            family = gaussian(), data = dd)
  expect_false(any(grepl("^cor_", variables(c0))))
  # the key written twice, and beside ||, is refused
  expect_error(frm(bf(y ~ x + (1 | b | gr(g1, id = "a"))),
                   family = gaussian(), data = dd),
               "names its correlation id twice", fixed = TRUE)
  expect_error(frm(bf(y ~ x + (1 + x || gr(g1, id = "a"))),
                   family = gaussian(), data = dd),
               "An |ID| key correlates", fixed = TRUE)
})

test_that("re_formula reads the term inside a covariance wrapper", {
  dd <- lab_data()
  fit <- frm(bf(y ~ x + diag(1 + x | g1) + (1 | g)), family = gaussian(),
             data = dd)
  nd <- dd[1:4, ]
  # diag(1 + x | g1) names the g1 term, as (1 + x | g1) does
  expect_identical(fitted(fit, newdata = nd, re_formula = ~ diag(1 + x | g1)),
                   fitted(fit, newdata = nd, re_formula = ~ (1 + x | g1)))
  expect_false(identical(
    fitted(fit, newdata = nd, re_formula = ~ diag(1 + x | g1)),
    fitted(fit, newdata = nd, re_formula = NA)))
})

test_that("a refused cor prior names a || term as it was written", {
  dd <- lab_data()
  dd$f <- factor(sample(c("a", "b", "c"), nrow(dd), TRUE))
  # the term's label is 1 + f | g, the correlated term it is not
  expect_error(frm(bf(y ~ x + (1 + f || g)), family = gaussian(), data = dd,
                   prior = set_prior("lkj(2)", class = "cor")),
               "no correlation parameter: (1 + f || g)", fixed = TRUE)
  expect_error(frm(bf(y ~ x + (1 + f | gr(g, cor = FALSE))),
                   family = gaussian(), data = dd,
                   prior = set_prior("lkj(2)", class = "cor")),
               "no correlation parameter: (1 + f | gr(g, cor = FALSE))",
               fixed = TRUE)
  # a term written diag() keeps its label and structure
  expect_error(frm(bf(y ~ x + diag(1 + f | g)), family = gaussian(),
                   data = dd, prior = set_prior("lkj(2)", class = "cor")),
               "no correlation parameter: 1 + f | g [diag]", fixed = TRUE)
})

test_that("variables() lists sd, then cor, then a smooth's sds", {
  # brms's reorder_pars(): sd, cor, df, sds, sdgp, lscale
  # (dev/prilabel-rev-log/order-brms.txt)
  dd <- lab_data()
  fit <- frm(bf(y ~ s(x) + (1 + x | g)), family = gaussian(), data = dd)
  v <- variables(fit)
  expect_identical(v[grepl("^(sd|cor|sds)_", v)],
                   c("sd_g__Intercept", "sd_g__x", "cor_g__Intercept__x",
                     "sds_sx_1"))
})

test_that("groups are listed in brms's order, sorted by name", {
  # brms's terms_re(): out[order(out$group), ] within a predictor
  # (dev/prilabel-rev-log/brmsfit.txt, fit A)
  dd <- lab_data()
  dd$k1 <- factor(sample(c("c:d", "e"), nrow(dd), TRUE))
  dd$k2 <- factor(sample(c("f", "h"), nrow(dd), TRUE))
  fit <- frm(bf(y ~ x + (1 + x | g1/g2) + (1 | g) + (1 | mm(m1, m2)) +
                  (1 | k1/k2)), family = gaussian(), data = dd)
  grp <- c("g", "g1", "g1:g2", "k1", "k1:k2", "mmm1m2")
  expect_identical(names(ranef(fit)), grp)
  expect_identical(names(ngrps(fit)), grp)
  expect_identical(setdiff(names(VarCorr(fit)), "residual__"), grp)
  v <- variables(fit)
  expect_identical(sub("__.*", "", v[grepl("^sd_", v)]),
                   paste0("sd_", c("g", "g1", "g1", "g1:g2", "g1:g2", "k1",
                                   "k1:k2", "mmm1m2")))
  # a part holding a ':' keeps it: brms pastes the parts with '_'
  expect_identical(rownames(ranef(fit)[["k1:k2"]]),
                   c("c:d_f", "c:d_h", "e_f", "e_h"))
  nd <- dd[1:2, ]
  nd$k1 <- factor("e")
  nd$k2 <- factor("zz")
  expect_error(fitted(fit, newdata = nd),
               "grouping factor `k1:k2`: e_zz", fixed = TRUE)
})

test_that("re_formula names a group as brms writes it", {
  dd <- lab_data()
  fit <- frm(bf(y ~ x + (1 + x | g1/g2)), family = gaussian(), data = dd)
  nd <- dd[1:4, ]
  expect_false(identical(
    fitted(fit, newdata = nd, re_formula = ~ (1 + x | g1:g2)),
    fitted(fit, newdata = nd, re_formula = NA)))
  # reformulas's order is not a group of this fit, and brms keeps no
  # term for it; frmtmb refuses it rather than drop it silently
  expect_error(fitted(fit, newdata = nd, re_formula = ~ (1 + x | g2:g1)),
               "matches no group-level term", fixed = TRUE)
})

test_that("a block with no level map still reads brms's levels", {
  # a block whose parts could not be read keeps reformulas's label and
  # reads each ':' as a separator; that path once failed on any block
  # with more than one level (found by the fix-5 mutant,
  # dev/prilabel-prefix-run.R)
  bk <- list(levels = c("p1:q1", "p1:q2"), group_name = "g1:g2")
  expect_identical(frmtmb:::brms_levels(bk), c("p1_q1", "p1_q2"))
  bk$brms_levels <- c("p1:q1" = "p1_q1", "p1:q2" = "p1_q2")
  expect_identical(frmtmb:::brms_levels(bk), c("p1_q1", "p1_q2"))
})
