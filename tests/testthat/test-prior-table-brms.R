# The prior table's rows, bounds and refusals against brms 2.23.0's
# default_prior() and validate_prior() on the same formulas (brms run
# with no Stan compile, the lane priortab record of the round of
# 2026-11-11). brms's answers, per case:
#
#   y ~ z + (1 | g): sd lb = 0 on the class row only; sigma lb = 0
#   student(): nu lb = 1, sigma lb = 0
#   zero_inflated_poisson(): zi lb = 0, ub = 1
#   gen_extreme_value(): xi has no bound
#   mvbind(y, y2) ~ x + (1 + x | g): one cor row per group, no resp
#   bf(y ~ x + (1 + x | g), sigma ~ (1 + x | g)): one cor row, no dpar
#   mvbind(Y_a, Y_b): resp Ya and Yb
#   table + set_prior() and set_prior() + table: a brmsprior of both
#   lkj on rescor with coef or group: refused (Lrescor_y__y2)
#   lkj on cor with dpar or nlpar: refused (L_mub, L_a)
#
# Through 0.81.0 every lb and ub of these rows was NA, `+` failed in R's
# Ops dispatch, and the four refused spellings were accepted.

ptb_data <- function() {
  set.seed(1)
  n <- 120
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(1:10, each = 12)))
  d$y <- 1 + 0.3 * d$z + stats::rnorm(10)[d$g] * 0.3 +
    stats::rnorm(n, sd = 0.3)
  d$y2 <- d$y + stats::rnorm(n)
  d$Y_a <- d$y
  d$Y_b <- d$y2
  d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
  d$cnt <- stats::rpois(n, 3)
  d
}

ptb_row <- function(tab, ...) {
  key <- list(...)
  tab <- as.data.frame(tab)
  hit <- rep(TRUE, nrow(tab))
  for (nm in names(key)) hit <- hit & tab[[nm]] == key[[nm]]
  tab[hit, , drop = FALSE]
}

test_that("sd and sigma rows carry brms's lb = 0, on the class row only", {
  d <- ptb_data()
  tab <- default_prior(y ~ z + (1 | g), data = d)
  sd_rows <- ptb_row(tab, class = "sd")
  expect_identical(sd_rows$lb[sd_rows$group == ""], 0)
  expect_true(all(is.na(sd_rows$lb[sd_rows$group != ""])))
  expect_identical(ptb_row(tab, class = "sigma")$lb, 0)
  expect_true(all(is.na(ptb_row(tab, class = "b")$lb)))
  expect_true(all(is.na(tab$ub)))
})

test_that("a distributional class row carries its interval, as brms's", {
  d <- ptb_data()
  st <- default_prior(y ~ z, data = d, family = student())
  expect_identical(ptb_row(st, class = "nu")$lb, 1)
  zi <- default_prior(cnt ~ z, data = d, family = zero_inflated_poisson())
  expect_identical(c(ptb_row(zi, class = "zi")$lb,
                     ptb_row(zi, class = "zi")$ub), c(0, 1))
  # the guard's absent case: brms prints no bound on xi, whose link
  # keeps it above -1, and none on a dpar given a formula
  gv <- default_prior(y ~ z, data = d, family = gen_extreme_value())
  expect_true(is.na(ptb_row(gv, class = "xi")$lb))
  sg <- default_prior(bf(y ~ z, sigma ~ z), data = d)
  expect_true(all(is.na(ptb_row(sg, dpar = "sigma")$lb)))
})

test_that("a table passed back moves no fit and leaves it an ML fit", {
  d <- ptb_data()
  cases <- list(list(y ~ z + (1 | g), gaussian()),
                list(y ~ z + (1 | g), student()),
                # an sds row carried lb = 0 before this round, and the
                # fit read it as a prior
                list(y ~ s(x, k = 5), gaussian()))
  # a Student-t on near-gaussian data leaves nu undetermined, on both
  w <- "Standard errors are not available"
  for (cs in cases) {
    a <- allow_warnings(frm(cs[[1]], data = d, family = cs[[2]]), w)
    b <- allow_warnings(
      frm(cs[[1]], data = d, family = cs[[2]],
          prior = default_prior(cs[[1]], data = d, family = cs[[2]])), w)
    expect_lte(abs(a$opt$objective - b$opt$objective),
               1e-10 * abs(a$opt$objective))
    expect_identical(names(a$opt$par), names(b$opt$par))
    expect_null(b$prior)
  }
  # the guard's absent case: a bound that boxes something is kept
  cbox <- frm(y ~ z + (1 | g), data = d,
              prior = set_prior("", class = "sd", lb = 0.05))
  expect_s3_class(cbox$prior, "frmtmb_priorlist")
  expect_identical(unname(resolve_prior_input(cbox, cbox$prior)$lower),
                   log(0.05))
  # and a bounds-only specification the user wrote is kept even where
  # it boxes nothing: brms reads it as (flat) on that slot, which
  # frmtmb.sample's frm_sample(fit) stacks over its defaults
  for (cl in c("sd", "sigma")) {
    own <- frm(y ~ z + (1 | g), data = d,
               prior = set_prior("", class = cl, lb = 0))
    expect_s3_class(own$prior, "frmtmb_priorlist")
    # a table row the user marked as theirs is kept too
    tab <- default_prior(y ~ z + (1 | g), data = d)
    tab$source[tab$class == cl & tab$group == ""] <- "user"
    expect_identical(nrow(as.data.frame(
      frm(y ~ z + (1 | g), data = d, prior = tab)$prior)), 1L)
  }
})

test_that("the sampling route adds a bound its default carries", {
  local_frmtmb_registries()
  frmtmb_register_prior_defaults(function(spec, frame) {
    set_prior("student_t(3, 0, 2.5)", class = "sigma", lb = 0.5)
  })
  d <- ptb_data()
  ps <- default_prior(y ~ z, data = d, route = "sample")
  expect_identical(ptb_row(ps, class = "sigma")$lb, 0.5)
  expect_identical(ptb_row(default_prior(y ~ z, data = d),
                           class = "sigma")$lb, 0)
  # no attribute leaks into the prior column
  expect_null(attributes(ps$prior))
})

test_that("a table and set_prior() combine with + in either order", {
  d <- ptb_data()
  p <- default_prior(y ~ z + x, data = d)
  sp <- set_prior("normal(0, 1)", class = "b", coef = "z")
  ab <- expect_no_warning(p + sp)
  ba <- expect_no_warning(sp + p)
  expect_s3_class(ab, "frmtmb_priorlist")
  expect_s3_class(ba, "frmtmb_priorlist")
  expect_identical(nrow(as.data.frame(ab)), nrow(p) + 1L)
  expect_identical(attr(ab, "route"), "fit")
  expect_s3_class(p + p, "frmtmb_priorlist")
  # the table's own (flat) row for coef z and the added one name one
  # slot, which brms's validate_prior() refuses too
  expect_error(validate_prior(ab, y ~ z + x, data = d),
               "Duplicated prior specifications", class = "frmtmb_error")
  pe <- p[p$coef != "z", ]
  fit <- frm(y ~ z + x, data = d, prior = pe + sp)
  expect_identical(nrow(as.data.frame(fit$prior)), 1L)
})

test_that("rescor takes no coef and no group, as brms requires", {
  d <- ptb_data()
  fm <- bf(mvbind(y, y2) ~ z) + set_rescor(TRUE)
  expect_error(validate_prior(set_prior("lkj(2)", class = "rescor",
                                        coef = "y__y2"), fm, data = d),
               "takes no coef", class = "frmtmb_error")
  expect_error(validate_prior(set_prior("lkj(2)", class = "rescor",
                                        group = "g"), fm, data = d),
               "takes no group", class = "frmtmb_error")
  tab <- validate_prior(set_prior("lkj(2)", class = "rescor"), fm, data = d)
  expect_identical(ptb_row(tab, class = "rescor")$prior, "lkj(2)")
})

test_that("cor takes no dpar or nlpar; its rows are keyed by group", {
  d <- ptb_data()
  expect_error(validate_prior(set_prior("lkj(2)", class = "cor",
                                        dpar = "mub"),
                              bf(f ~ x + (1 | p | g)), data = d,
                              family = categorical()),
               "takes no dpar", class = "frmtmb_error")
  fs <- bf(y ~ x + (1 + x | g), sigma ~ (1 + x | g))
  expect_error(validate_prior(set_prior("lkj(2)", class = "cor",
                                        dpar = "sigma"), fs, data = d),
               "takes no dpar", class = "frmtmb_error")
  fnl <- bf(y ~ a + b * x, a ~ 1 + (1 + x | g), b ~ 1, nl = TRUE)
  expect_error(validate_prior(set_prior("lkj(2)", class = "cor",
                                        nlpar = "a"), fnl, data = d),
               "takes no nlpar", class = "frmtmb_error")
  for (f in list(fs, fnl, bf(mvbind(y, y2) ~ x + (1 + x | g)))) {
    cr <- ptb_row(default_prior(f, data = d), class = "cor")
    expect_identical(cr$group, c("", "g"))
    expect_true(all(cr$dpar == "" & cr$nlpar == "" & cr$resp == ""))
  }
  # the group spelling reaches both blocks, as brms's cor g does
  tab <- validate_prior(set_prior("lkj(4)", class = "cor", group = "g"), fs,
                        data = d)
  expect_identical(ptb_row(tab, class = "cor", group = "g")$prior, "lkj(4)")
  # resp stays frmtmb's own spelling (the decision of 2026-09-24)
  mv <- bf(mvbind(y, y2) ~ x + (1 + x | g))
  tv <- validate_prior(set_prior("lkj(5)", class = "cor", group = "g",
                                 resp = "y2"), mv, data = d)
  expect_identical(ptb_row(tv, class = "cor", resp = "y2")$prior, "lkj(5)")
})

test_that("a response is listed by its Stan name, as in brms", {
  d <- ptb_data()
  f <- bf(mvbind(Y_a, Y_b) ~ z) + set_rescor(FALSE)
  tab <- default_prior(f, data = d)
  expect_setequal(unique(tab$resp), c("Ya", "Yb"))
  for (r in c("Y_a", "Ya")) {
    tv <- validate_prior(set_prior("normal(0, 1)", class = "b", resp = r),
                         f, data = d)
    expect_identical(ptb_row(tv, class = "b", coef = "", resp = "Ya")$prior,
                     "normal(0, 1)")
    expect_false(any(tv$resp == "Y_a"))
  }
  # the two names are one slot: two specifications are refused, where
  # one of them used to be dropped in silence
  expect_error(frm(f, data = d,
                   prior = set_prior("normal(0, 1)", class = "b",
                                     resp = "Y_a") +
                     set_prior("normal(5, 0.01)", class = "b",
                               resp = "Ya")),
               "Duplicated prior specifications", class = "frmtmb_error")
  # the table passed back reaches the response it names
  tab$prior[tab$class == "b" & tab$coef == "z" & tab$resp == "Yb"] <-
    "normal(0, 0.1)"
  fit <- frm(f, data = d, prior = tab)
  pr <- as.data.frame(fit$prior)
  expect_identical(pr$resp[pr$prior == "normal(0, 0.1)"], "Yb")
})

test_that("prior = refusals say what was given", {
  d <- ptb_data()
  sp <- set_prior("normal(0, 1)")
  expect_error(frm(y ~ z, data = d, prior = list(sp)),
               "Pass the specification itself", class = "frmtmb_error")
  expect_error(frm(y ~ z, data = d, prior = list(sp, sp)),
               "sp1 + sp2", fixed = TRUE, class = "frmtmb_error")
  expect_error(frm(y ~ z, data = d, prior = list()), "an empty list",
               class = "frmtmb_error")
  expect_error(frm(y ~ z, data = d, prior = "flat"),
               "frm_sample()'s spelling", fixed = TRUE,
               class = "frmtmb_error")
  expect_error(frm(y ~ z, data = d, prior = "normal(0, 1)"),
               "set_prior(\"normal(0, 1)\", class = \"b\")", fixed = TRUE,
               class = "frmtmb_error")
})

test_that("a table cut to fewer columns is refused for what it lost", {
  d <- ptb_data()
  pe <- default_prior(y ~ z + x, data = d)
  pe$prior[pe$class == "b" & pe$coef == "z"] <- "normal(0, 1)"
  expect_error(frm(y ~ z + x, data = d, prior = pe[, "prior", drop = FALSE]),
               "without the columns `class`", class = "frmtmb_error")
  expect_error(frm(y ~ z + x, data = d, prior = pe[, c("prior", "class")]),
               "without the columns `coef`", class = "frmtmb_error")
  # rows, or the bound and source columns, can go
  keep <- c("prior", "class", "coef", "group", "resp", "dpar", "nlpar")
  a <- frm(y ~ z + x, data = d, prior = pe[pe$class == "b", ])
  b <- frm(y ~ z + x, data = d, prior = pe[, keep])
  expect_identical(a$opt$par, b$opt$par)
  # a plain data frame keeps set_prior()'s defaults, as.brmsprior()'s
  # documented reading
  expect_s3_class(as.brmsprior(data.frame(prior = "normal(0, 1)",
                                          coef = "z")), "frmtmb_priorlist")
})

test_that("a uniform on a slope with no bounds boxes it to its support", {
  d <- ptb_data()
  d$y3 <- d$y + 0.5 * d$z
  w <- c("no natural lower bound", "no natural upper bound",
         "Standard errors are not available")
  a <- allow_warnings(frm(y3 ~ z, data = d,
                          prior = set_prior("uniform(0, 0.2)", class = "b",
                                            coef = "z")),
                      w, require = w[1:2])
  b <- allow_warnings(frm(y3 ~ z, data = d,
                          prior = set_prior("uniform(0, 0.2)", class = "b",
                                            coef = "z", lb = 0, ub = 0.2)),
                      w)
  zof <- function(f) unname(fixef(f)["z", "Estimate"])
  expect_identical(a$opt$convergence, 0L)
  expect_identical(a$opt$par, b$opt$par)
  expect_identical(zof(a), 0.2)
  # an interior mode is the unpenalized one: a uniform is a constant
  # inside its support
  ml <- frm(y3 ~ z, data = d)
  ai <- allow_warnings(frm(y3 ~ z, data = d,
                           prior = set_prior("uniform(0, 2)", class = "b",
                                             coef = "z")), w)
  expect_lte(abs(zof(ai) - zof(ml)), 1e-5 * abs(zof(ml)))
  # the guard's absent case: a density on the whole line boxes nothing
  ri <- resolve_prior_input(ml, set_prior("normal(0, 1)", class = "b",
                                          coef = "z"))
  expect_length(ri$lower, 0L)
  expect_length(ri$upper, 0L)
  ru <- resolve_prior_input(ml, set_prior("uniform(0, 2)", class = "b",
                                          coef = "z"))
  expect_identical(unname(c(ru$lower[["z"]], ru$upper[["z"]])), c(0, 2))
})

test_that("the box reaches a dpar's own class and an uncentered intercept", {
  d <- ptb_data()
  w <- c("no natural lower bound", "no natural upper bound",
         "Standard errors are not available")
  pair <- function(f, sp, lb, ub) {
    a <- allow_warnings(frm(f, data = d, prior = sp), w)
    s2 <- sp
    s2$lb <- lb
    s2$ub <- ub
    b <- allow_warnings(frm(f, data = d, prior = s2), w)
    expect_identical(a$opt$convergence, 0L)
    expect_identical(a$opt$par, b$opt$par)
    a
  }
  a <- pair(y ~ z, set_prior("uniform(0.6, 1)", class = "sigma"), 0.6, 1)
  expect_identical(unname(sigma(a)), 0.6)
  a <- pair(bf(y ~ 0 + Intercept + z),
            set_prior("uniform(0, 0.5)", class = "b", coef = "Intercept"),
            0, 0.5)
  expect_identical(unname(fixef(a)["Intercept", "Estimate"]), 0.5)
  # an exponential is positive at 0, so it is boxed there as well
  d$yneg <- d$y - 0.8 * d$z
  ex <- allow_warnings(frm(yneg ~ z, data = d,
                           prior = set_prior("exponential(1)", class = "b",
                                             coef = "z")), w)
  expect_identical(ex$opt$convergence, 0L)
  expect_true(is.finite(ex$opt$objective))
  expect_identical(unname(fixef(ex)["z", "Estimate"]), 0)
})

test_that("zero-density ends get no box; an Inf stop is not converged", {
  # the review's case: beta(2, 2) on a slope the data put below 0. Its
  # mode is inside (0, 1), near 0.0068, where lb = 1e-8 finds it; a
  # closed box at 0 stopped there with the objective Inf and code 0
  d <- ptb_data()
  d$yneg <- d$y - 0.8 * d$z
  w <- c("no natural lower bound", "no natural upper bound")
  sp <- function(...) set_prior("beta(2, 2)", class = "b", coef = "z", ...)
  zof <- function(f) unname(fixef(f)["z", "Estimate"])
  m8 <- frm(yneg ~ z, data = d, prior = sp(lb = 1e-8, ub = 1))
  m10 <- frm(yneg ~ z, data = d, prior = sp(lb = 1e-10, ub = 1))
  expect_identical(m8$opt$convergence, 0L)
  expect_true(zof(m8) > 1e3 * 1e-8)
  expect_lte(abs(m8$opt$objective - m10$opt$objective),
             1e-8 * abs(m8$opt$objective))
  # with no bound it is not boxed: the fit fails as before, loudly,
  # rather than stop at the end
  expect_error(allow_warnings(frm(yneg ~ z, data = d, prior = sp()),
                              c(w, "NA/NaN function evaluation")),
               "optimizer failed", class = "frmtmb_error")
  # and a mode inside is found with no bound, as before
  pos <- allow_warnings(frm(y ~ z, data = d, prior = sp()), w)
  expect_identical(pos$opt$convergence, 0L)
  expect_true(zof(pos) > 0 && zof(pos) < 1)
  # the user's closed box at the end: a stop at an infinite objective
  # is reported as not converged
  at0 <- allow_warnings(frm(yneg ~ z, data = d,
                            prior = sp(lb = 0, ub = 1)),
                        c("the objective is Inf", "gradient at the reported"),
                        require = "the objective is Inf")
  expect_identical(at0$opt$convergence, 1L)
  expect_false(is.finite(at0$opt$objective))
})

test_that("a link whose inverse decays to 0 gives lb = 0, as brms", {
  d <- ptb_data()
  for (lk in c("softplus", "squareplus")) {
    tab <- as.data.frame(default_prior(
      y ~ z, data = d, family = brmsfamily("gaussian", link_sigma = lk)))
    expect_identical(tab$lb[tab$class == "sigma"], 0)
  }
  # the guard's absent case: an identity link leaves sigma unbounded
  tab <- as.data.frame(default_prior(
    y ~ z, data = d, family = brmsfamily("gaussian",
                                         link_sigma = "identity")))
  expect_true(is.na(tab$lb[tab$class == "sigma"]))
})
