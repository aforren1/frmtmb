# A class "sd" prior narrowed to one coefficient, as brms 2.23.0 lists
# and honors it: default_prior() has a row per coefficient below each
# group row, and set_prior(class = "sd", group = "g", coef = "x") puts
# the density on sd_g__x alone. The expected tables are brms's own
# (dev/prilabel-log/brms-ref.txt, from dev/prilabel-brms-ref.R); the
# gated test at the end compares against brms directly.

sdcoef_data <- function(n = 240, seed = 11) {
  set.seed(seed)
  g <- factor(sample(sprintf("g%02d", 1:20), n, TRUE))
  h <- factor(sample(sprintf("h%d", 1:8), n, TRUE))
  f <- factor(sample(c("a", "b", "c"), n, TRUE))
  x <- rnorm(n)
  u0 <- rnorm(20, 0, 0.8)[as.integer(g)]
  u1 <- rnorm(20, 0, 0.4)[as.integer(g)]
  y <- 1 + 0.5 * x + u0 + u1 * x + rnorm(8, 0, 0.5)[as.integer(h)] +
    rnorm(n, 0, 0.7)
  data.frame(y, x, f, g, h)
}

sd_rows <- function(form, data) {
  p <- as.data.frame(default_prior(form, data = data, family = gaussian()))
  p <- p[p$class == "sd", c("coef", "group", "dpar")]
  rownames(p) <- NULL
  p
}

test_that("default_prior() lists an sd row per coefficient, as brms", {
  dd <- sdcoef_data()
  expect_identical(
    sd_rows(bf(y ~ x + (1 + x | g)), dd),
    data.frame(coef = c("", "", "Intercept", "x"),
               group = c("", "g", "g", "g"), dpar = ""))
  # brms sorts the coefficient rows with order(), so the collation
  # decides: fb, fc, Intercept in en_US, Intercept first in C
  expect_identical(
    sd_rows(bf(y ~ x + (f | g)), dd),
    data.frame(coef = c("", "", sort(c("fb", "fc", "Intercept"))),
               group = c("", "g", "g", "g", "g"), dpar = ""))
  # a block on sigma: its rows follow the location's, each prefix with
  # its own class-wide row
  expect_identical(
    sd_rows(bf(y ~ x + (1 | g), sigma ~ (1 + x | g)), dd),
    data.frame(coef = c("", "", "Intercept", "", "", "Intercept", "x"),
               group = c("", "g", "g", "", "g", "g", "g"),
               dpar = c("", "", "", "sigma", "sigma", "sigma", "sigma")))
  # || is brms's gr(g, cor = FALSE), with the same rows as |
  expect_identical(sd_rows(bf(y ~ x + (1 + x || g)), dd),
                   sd_rows(bf(y ~ x + (1 + x | g)), dd))
})

test_that("a per-coefficient sd prior reaches that sd alone", {
  dd <- sdcoef_data()
  form <- bf(y ~ x + (1 + x | g) + (1 | h))
  pr <- set_prior("normal(0, 0.5)", class = "sd", group = "g", coef = "x")
  f0 <- frm(form, family = gaussian(), data = dd)
  f1 <- frm(form, family = gaussian(), data = dd, prior = pr)
  bk <- Filter(function(b) identical(b$group_name, "g"),
               f1$frame$re_blocks)[[1L]]
  th <- names(f1$obj$par) == "theta"
  pos <- which(th)
  # theta of the g block: log sd(Intercept), log sd(x), then the
  # correlation; h's log sd follows
  i_int <- pos[bk$theta_idx[1L]]
  i_x <- pos[bk$theta_idx[2L]]
  i_cor <- pos[bk$theta_idx[3L]]
  p <- f1$obj$par
  pen <- function(q) as.numeric(f1$obj$fn(q) - f0$obj$fn(q))
  shift <- function(q, i, by) {
    q[i] <- q[i] + by
    q
  }
  # the penalty is the negative log density of sd(x) with the
  # log-Jacobian of theta = log sd; constants cancel in a difference
  ld <- function(t) stats::dnorm(exp(t), 0, 0.5, log = TRUE) + t
  for (by in c(-0.4, 0.3)) {
    expect_equal(pen(shift(p, i_x, by)) - pen(p),
                 -(ld(p[[i_x]] + by) - ld(p[[i_x]])))
    # the intercept's sd, the correlation and the other group carry
    # none of it
    for (i in c(i_int, i_cor, pos[length(pos)])) {
      expect_equal(pen(shift(p, i, by)), pen(p))
    }
  }
  # and the prior table says the same, row by row
  tab <- as.data.frame(prior_summary(f1))
  sdx <- tab$class == "sd" & tab$group == "g" & tab$coef == "x"
  expect_identical(tab$prior[sdx], "normal(0, 0.5)")
  other <- tab$class == "sd" & tab$coef == "Intercept"
  expect_true(all(tab$prior[other] == "(flat)"))
})

test_that("a coef prior wins over its group's whatever the order", {
  dd <- sdcoef_data()
  form <- bf(y ~ x + (1 + x | g))
  grp <- set_prior("normal(0, 2)", class = "sd", group = "g")
  cf <- set_prior("normal(0, 0.5)", class = "sd", group = "g", coef = "x")
  for (pr in list(grp + cf, cf + grp)) {
    des <- frmtmb:::prior_design(form, dd, gaussian(), list())
    ent <- frmtmb:::resolve_priorlist(des, pr)$entries
    bk <- des$frame$re_blocks[[1L]]
    sc <- vapply(ent, function(e) e$dist$scale, 0)
    at <- as.integer(vapply(ent, function(e) e$idx, 0))
    # brms 2.23.0: normal(0, 2) on sd_1[1] and normal(0, 0.5) on sd_1[2]
    expect_identical(sc[match(bk$theta_idx[1:2], at)], c(2, 0.5))
  }
})

test_that("per-coefficient sd priors cover brms's spellings of a block", {
  dd <- sdcoef_data()
  one <- function(form, ...) {
    des <- frmtmb:::prior_design(form, dd, gaussian(), list())
    pr <- set_prior("normal(0, 0.5)", class = "sd", ...)
    ent <- frmtmb:::resolve_priorlist(des, pr)$entries
    as.integer(vapply(ent, function(e) e$idx, 0))
  }
  # a factor's own column, a || term, gr(cor = FALSE), a dpar block
  expect_identical(one(bf(y ~ x + (f | g)), group = "g", coef = "fb"), 2L)
  expect_identical(one(bf(y ~ x + (1 + x || g)), group = "g", coef = "x"),
                   2L)
  expect_identical(one(bf(y ~ x + (1 + x | gr(g, cor = FALSE))),
                       group = "g", coef = "x"), 2L)
  # (1 | g) on mu is theta 1; sigma's block is theta 2 to 4
  expect_identical(one(bf(y ~ x + (1 | g), sigma ~ (1 + x | g)),
                       group = "g", coef = "x", dpar = "sigma"), 3L)
  expect_identical(one(bf(y ~ x + (1 | g), sigma ~ (1 + x | g)),
                       group = "g", coef = "Intercept"), 1L)
  # brms puts a coef with no dpar on the location's column alone, even
  # in a block that spans mu and sigma (dev/prilabel-log/brms-ref3.txt)
  expect_identical(one(bf(y ~ x + (1 + x | p | g), sigma ~ (1 | p | g)),
                       group = "g", coef = "Intercept"), 1L)
})

test_that("a per-coefficient sd prior that names nothing is refused", {
  dd <- sdcoef_data()
  form <- bf(y ~ x + (1 + x | g))
  # brms: "The following priors do not correspond to any model
  # parameter: sd_x"
  expect_error(frm(form, family = gaussian(), data = dd,
                   prior = set_prior("normal(0, 1)", class = "sd",
                                     coef = "x")),
               "needs a `group`", fixed = TRUE)
  # brms: "... sd_g__z"
  expect_error(frm(form, family = gaussian(), data = dd,
                   prior = set_prior("normal(0, 1)", class = "sd",
                                     group = "g", coef = "z")),
               "group \"g\": \"Intercept\", \"x\"", fixed = TRUE)
  # a structure whose columns share one sd has none per coefficient
  expect_error(frm(bf(y ~ x + homdiag(1 + x | g)), family = gaussian(),
                   data = dd,
                   prior = set_prior("normal(0, 1)", class = "sd",
                                     group = "g", coef = "x")),
               "one shared standard deviation", fixed = TRUE)
  # and where each column has its own the same prior is honored
  for (form2 in list(bf(y ~ x + diag(1 + x | g)),
                     bf(y ~ x + cs(1 + x | g)))) {
    des <- frmtmb:::prior_design(form2, dd, gaussian(), list())
    pr <- set_prior("normal(0, 1)", class = "sd", group = "g", coef = "x")
    ent <- frmtmb:::resolve_priorlist(des, pr)$entries
    expect_identical(as.integer(vapply(ent, function(e) e$idx, 0)), 2L)
  }
  # the shared sd itself still takes a group prior
  des <- frmtmb:::prior_design(bf(y ~ x + homdiag(1 + x | g)), dd,
                               gaussian(), list())
  ent <- frmtmb:::resolve_priorlist(
    des, set_prior("normal(0, 1)", class = "sd", group = "g"))$entries
  expect_identical(as.integer(vapply(ent, function(e) e$idx, 0)), 1L)
})

test_that("brms's per-coefficient sd row translates and lands", {
  skip_if_not_installed("brms")
  dd <- sdcoef_data()
  pl <- frmtmb:::as_priorlist(brms::prior(exponential(1), class = "sd",
                                          group = "g", coef = "x"))
  des <- frmtmb:::prior_design(bf(y ~ x + (1 + x | g)), dd, gaussian(),
                               list())
  ent <- frmtmb:::resolve_priorlist(des, pl)$entries
  expect_identical(as.integer(vapply(ent, function(e) e$idx, 0)), 2L)
  # every sd row of brms's table is one of frmtmb's, in brms's order
  for (form in list(y ~ x + (1 + x | g) + (1 | h), y ~ x + (f | g))) {
    b <- as.data.frame(brms::default_prior(brms::bf(form), data = dd))
    b <- b[b$class == "sd", c("coef", "group", "dpar")]
    rownames(b) <- NULL
    expect_identical(sd_rows(bf(form), dd), b)
  }
})

test_that("a shared-sd block refuses a coef of its own prefix only", {
  # lane prilabel's review, finding 2 (dev/prilabel-rev-probe.R): brms
  # puts a coef prior with no dpar on the location's block alone, so a
  # homdiag() block on sigma is no reason to refuse it
  dd <- sdcoef_data()
  pr <- set_prior("normal(0, 0.5)", class = "sd", group = "g", coef = "x")
  des <- frmtmb:::prior_design(bf(y ~ x + (1 + x | g),
                                  sigma ~ homdiag(1 + x | g)),
                               dd, gaussian(), list())
  ent <- frmtmb:::resolve_priorlist(des, pr)$entries
  # mu's block is theta 1 to 3 (two sds, one correlation): x is theta 2
  expect_identical(as.integer(vapply(ent, function(e) e$idx, 0)), 2L)
  # the same coef on sigma's own block names its shared sd
  expect_error(frmtmb:::resolve_priorlist(
    des, set_prior("normal(0, 0.5)", class = "sd", group = "g", coef = "x",
                   dpar = "sigma")),
    "one shared standard deviation", fixed = TRUE)
  # and the other way round: the shared block on mu, a coef on sigma's
  des2 <- frmtmb:::prior_design(bf(y ~ x + homdiag(1 + x | g),
                                   sigma ~ (1 + x | g)),
                                dd, gaussian(), list())
  ent2 <- frmtmb:::resolve_priorlist(
    des2, set_prior("normal(0, 0.5)", class = "sd", group = "g",
                    coef = "x", dpar = "sigma"))$entries
  # mu's homdiag block is theta 1; sigma's x is theta 3
  expect_identical(as.integer(vapply(ent2, function(e) e$idx, 0)), 3L)
  expect_error(frmtmb:::resolve_priorlist(des2, pr),
               "one shared standard deviation", fixed = TRUE)
})

test_that("validate_prior() lists the sd rows in brms's order", {
  # brms's validate_prior(): order(class, group, resp, dpar, nlpar, coef),
  # every class-wide sd row first, then group by group, the responses
  # inside each (lane prilabel's review, finding 8;
  # dev/prilabel-rev-log/brmsfit.txt, fit B)
  dd <- sdcoef_data()
  dd$y2 <- dd$y + stats::rnorm(nrow(dd))
  vp <- as.data.frame(validate_prior(
    set_prior("normal(0, 0.5)", class = "sd", group = "g", coef = "x",
              resp = "y"),
    mvbf(bf(y ~ x + (1 + x | g)), bf(y2 ~ x + (1 | g))), data = dd,
    family = gaussian()))
  sd <- vp[vp$class == "sd", c("coef", "group", "resp")]
  rownames(sd) <- NULL
  expect_identical(sd, data.frame(
    coef = c("", "", "", "Intercept", "x", "", "Intercept"),
    group = c("", "", "g", "g", "g", "g", "g"),
    resp = c("y", "y2", "y", "y", "y", "y2", "y2")))
  expect_identical(vp$prior[vp$class == "sd" & vp$coef == "x"],
                   "normal(0, 0.5)")
})
