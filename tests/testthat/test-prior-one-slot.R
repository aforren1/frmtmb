# One prior slot per random-effect and residual-correlation parameter,
# as in brms. brms 2.23.0 has no class "theta": validate_prior() refuses
# `set_prior(class = "theta", coef = "theta_4")` ("do not correspond to
# any model parameter"), so a standard deviation is reached by class
# "sd", a correlation by class "cor" (lkj() on the whole matrix, and no
# prior on one correlation alone), a residual correlation by "rescor",
# and a gp() term's parameters by "sdgp" and "lscale". brms accepts
# lb/ub on class "cor" and then writes a Stan program that does not
# parse, so it has no bound on a correlation either.
#
# Through 0.76.0 class "theta" was a second spelling of each, and the
# one written later replaced the other with no message (seed 1, block
# (1 + x + z | g), 8 groups of 10): `cor lkj(2)` then
# `theta theta_4 normal(0, 1)` kept only the normal on position 4 and
# dropped the LKJ on 5 and 6, so a MAP fit took theta_5 and theta_6 to
# 558.85 and 1740.96; in frm_sample() one theta_4 prior retired the
# default lkj(1) on all three correlations; `sd` with coef = "x" plus
# `theta theta_2` dropped the sd prior. Class "theta" keeps the
# parameters no brms class reaches: a toep() correlation and a
# registered structure's own parameter.

ros_data <- function() {
  set.seed(1)
  G <- 8L
  n <- 10L
  d <- data.frame(g = factor(rep(seq_len(G), each = n)),
                  x = stats::rnorm(G * n), z = stats::rnorm(G * n),
                  t = rep(seq_len(n), G))
  d$tf <- factor(d$t)
  d$y <- 1 + 0.5 * d$x + stats::rnorm(G, 0, 0.7)[d$g] +
    stats::rnorm(G, 0, 0.4)[d$g] * d$x + stats::rnorm(G * n, 0, 0.6)
  d$y2 <- 0.5 * d$y + stats::rnorm(G * n)
  d
}

# the internal parameters a prior puts a density or a bound on, by name
ros_reached <- function(formula, data, pl, family = NULL) {
  des <- frmtmb:::prior_design(formula, data, family, list())
  r <- frmtmb:::resolve_priorlist(des, pl)
  pt <- des$frame[["par_template"]]
  nm <- unlist(lapply(r$entries, function(e) {
    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
  }))
  sort(unique(c(nm, names(r$lower), names(r$upper))))
}

ros_msg <- function(expr) {
  tryCatch({
    force(expr)
    ""
  }, error = conditionMessage)
}

test_that("class theta on a random-effect correlation is refused", {
  d <- ros_data()
  f <- bf(y ~ x + z + (1 + x + z | g))
  lkj <- set_prior("lkj(2)", class = "cor")
  th4 <- set_prior("normal(0, 1)", class = "theta", coef = "theta_4")
  # in either order: one parameter, one slot
  for (pl in list(lkj + th4, th4 + lkj, th4)) {
    m <- ros_msg(validate_prior(pl, f, data = d))
    expect_match(m, "internal parameter of a random-effect correlation",
                 fixed = TRUE)
    expect_match(m, "set_prior(..., class = \"cor\", group = \"g\")",
                 fixed = TRUE)
  }
  # a box on one correlation, which brms has no prior for either
  expect_match(ros_msg(validate_prior(
    set_prior("", class = "theta", coef = "theta_5", lb = -1, ub = 1), f,
    data = d)), "no prior on one correlation alone", fixed = TRUE)
  # the class spelling still reaches the whole matrix
  expect_identical(ros_reached(f, d, lkj), c("theta_4", "theta_5",
                                             "theta_6"))
})

test_that("class theta on a standard deviation names class sd", {
  d <- ros_data()
  f <- bf(y ~ x + z + (1 + x + z | g))
  sdx <- set_prior("student_t(3, 0, 1)", class = "sd", group = "g",
                   coef = "x")
  m <- ros_msg(validate_prior(
    sdx + set_prior("normal(0, 1)", class = "theta", coef = "theta_2"), f,
    data = d))
  expect_match(m, paste0("set_prior(..., class = \"sd\", group = \"g\", ",
                         "coef = \"x\")"), fixed = TRUE)
  expect_match(m, "on the standard deviation itself", fixed = TRUE)
  # the bare position, frmtmb's older spelling, is the same parameter
  expect_match(ros_msg(validate_prior(
    set_prior("", class = "theta", coef = "2", lb = -2), f, data = d)),
    "class = \"sd\"", fixed = TRUE)
  expect_identical(ros_reached(f, d, sdx), "theta_2")
  # and class theta with no coef has no parameter left here
  expect_match(ros_msg(validate_prior(
    set_prior("normal(0, 2)", class = "theta"), f, data = d)),
    "every random-effect covariance parameter of this model has a class",
    fixed = TRUE)
})

test_that("a MAP fit takes the class priors and nothing else", {
  # through 0.76.0 the theta row replaced the LKJ on positions 5 and 6,
  # and this fit ran them off to 558.85 and 1740.96
  d <- ros_data()
  f <- bf(y ~ x + z + (1 + x + z | g))
  expect_error(
    frm(f, data = d, family = gaussian(),
        prior = set_prior("lkj(2)", class = "cor") +
          set_prior("normal(0, 1)", class = "theta", coef = "theta_4")),
    "has a class of its own", class = "frmtmb_error")
})

test_that("class theta on rescor, sdgp and a smooth's sd is refused", {
  d <- ros_data()
  fr <- bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE)
  m <- ros_msg(validate_prior(
    set_prior("lkj(2)", class = "rescor") +
      set_prior("normal(0, 1)", class = "theta", coef = "thetar_1"),
    fr, data = d))
  expect_match(m, "set_prior(..., class = \"rescor\")", fixed = TRUE)
  expect_match(ros_msg(validate_prior(
    set_prior("", class = "theta", coef = "thetar_1", lb = 0), fr,
    data = d)), "internal parameter of a residual correlation",
    fixed = TRUE)
  expect_identical(ros_reached(fr, d, set_prior("lkj(2)", class = "rescor")),
                   "thetar_1")

  m <- ros_msg(validate_prior(
    set_prior("normal(0, 1)", class = "theta", coef = "theta_1"),
    bf(y ~ gp(x)), data = d))
  expect_match(m, "set_prior(..., class = \"sdgp\", coef = \"gpx\")",
               fixed = TRUE)
  m <- ros_msg(validate_prior(
    set_prior("normal(0, 1)", class = "theta", coef = "theta_2"),
    bf(y ~ gp(x)), data = d))
  expect_match(m, "set_prior(..., class = \"lscale\", coef = \"gpx\")",
               fixed = TRUE)
  m <- ros_msg(validate_prior(
    set_prior("normal(0, 1)", class = "theta", coef = "theta_1"),
    bf(y ~ s(x)), data = d))
  expect_match(m, "set_prior(..., class = \"sd\", group = \"s(x)\")",
               fixed = TRUE)
})

test_that("a correlation of a structure brms lacks follows class cor", {
  d <- ros_data()
  # ar1() holds one shared sd and one correlation, which lkj() reaches
  f <- bf(y ~ x + ar1(0 + tf | g))
  expect_match(ros_msg(validate_prior(
    set_prior("normal(0, 1)", class = "theta", coef = "theta_2"), f,
    data = d)), "class = \"cor\", group = \"g\"", fixed = TRUE)
  expect_identical(ros_reached(f, d, set_prior("lkj(2)", class = "cor")),
                   "theta_2")
})

test_that("class theta keeps the parameters no class reaches", {
  d <- ros_data()
  # toep(): one sd per column (theta_1 to theta_10) and banded
  # correlations no LKJ density fits (theta_11 to theta_19)
  f <- bf(y ~ x + toep(0 + tf | g))
  expect_identical(ros_reached(f, d, set_prior("normal(0, 1)",
                                               class = "theta",
                                               coef = "theta_11")),
                   "theta_11")
  expect_identical(ros_reached(f, d, set_prior("", class = "theta",
                                               coef = "theta_12",
                                               lb = -1, ub = 1)),
                   "theta_12")
  # the class-wide spelling covers exactly those
  expect_identical(ros_reached(f, d, set_prior("normal(0, 2)",
                                               class = "theta")),
                   sort(paste0("theta_", 11:19)))
  expect_match(ros_msg(validate_prior(
    set_prior("normal(0, 1)", class = "theta", coef = "theta_3"), f,
    data = d)), "class = \"sd\", group = \"g\", coef = \"tf3\"",
    fixed = TRUE)
  tab <- default_prior(f, data = d)
  th <- tab[tab$class == "theta", ]
  expect_identical(th$coef, c("", paste0("theta_", 11:19)))
})

test_that("a registered structure's own parameter keeps class theta", {
  local_frmtmb_registries()
  # xrho(0 + t | g): one standard deviation (sd_idx) and an exchangeable
  # correlation tanh(theta_2) that no brms class names
  frmtmb_register_covstruct(
    "xrho",
    lpdf = function(b, theta, blk) {
      d <- blk[["dim"]]
      r <- tanh(theta[2])
      S <- exp(2 * theta[1]) * ((1 - r) * diag(d) + r)
      dim(b) <- c(d, length(b) %/% d)
      sum(RTMB::dmvnorm(t(b), 0, S, log = TRUE))
    },
    vcov = function(theta, blk) {
      d <- blk[["dim"]]
      r <- tanh(theta[2])
      V <- exp(2 * theta[1]) * ((1 - r) * diag(d) + r)
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = 2L, sd_idx = 1L)
  d <- ros_data()
  d$t2 <- factor(d$t %% 2)
  f <- bf(y ~ x + xrho(0 + t2 | g))
  expect_identical(ros_reached(f, d, set_prior("normal(0, 1)",
                                               class = "theta",
                                               coef = "theta_2")),
                   "theta_2")
  expect_identical(ros_reached(f, d, set_prior("normal(0, 1)",
                                               class = "theta")),
                   "theta_2")
  expect_match(ros_msg(validate_prior(
    set_prior("normal(0, 1)", class = "theta", coef = "theta_1"), f,
    data = d)), "class = \"sd\", group = \"g\"", fixed = TRUE)
  tab <- default_prior(f, data = d)
  expect_identical(tab$coef[tab$class == "theta"], c("", "theta_2"))
})

test_that("default_prior() lists no class theta row a class reaches", {
  d <- ros_data()
  tab <- default_prior(bf(y ~ x + z + (1 + x + z | g)), data = d)
  expect_false(any(tab$class == "theta"))
  tab <- default_prior(bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), data = d)
  expect_false(any(tab$class == "theta"))
  # brms 2.23.0's own rows for this model; the order is the session's
  # collation, which testthat sets to C
  tab <- default_prior(bf(y ~ x + (1 + x | g)), data = d)
  expect_setequal(paste(tab$class, tab$coef, tab$group),
                  c("b  ", "b x ", "cor  ", "cor  g", "Intercept  ",
                    "sd  ", "sd  g", "sd Intercept g", "sd x g",
                    "sigma  "))
})

test_that("a centered intercept has one slot, class Intercept", {
  d <- ros_data()
  m <- ros_msg(validate_prior(set_prior("normal(0, 1)", class = "b",
                                        coef = "Intercept"),
                              bf(y ~ x), data = d))
  expect_match(m, "has one prior slot, class = \"Intercept\"", fixed = TRUE)
  m <- ros_msg(validate_prior(set_prior("normal(0, 1)", class = "b",
                                        coef = "Intercept", dpar = "sigma"),
                              bf(y ~ x, sigma ~ x), data = d))
  expect_match(m, "class = \"Intercept\", dpar = \"sigma\"", fixed = TRUE)
  # an intercept that is a class b coefficient keeps that spelling
  expect_identical(ros_reached(bf(y ~ 0 + Intercept + x), d,
                               set_prior("normal(0, 1)", class = "b",
                                         coef = "Intercept")),
                   "(Intercept)")
  expect_identical(ros_reached(bf(y ~ x, center = FALSE), d,
                               set_prior("normal(0, 1)", class = "b",
                                         coef = "Intercept")),
                   "(Intercept)")
  expect_identical(ros_reached(bf(y ~ x), d,
                               set_prior("normal(0, 1)",
                                         class = "Intercept")),
                   "(Intercept)")
})

test_that("Intercept with nlpar is the slot of class b, coef Intercept", {
  d <- ros_data()
  d$yn <- 2 * exp(0.3 * d$x) + stats::rnorm(nrow(d), 0, 0.3)
  bi <- set_prior("normal(9, 1)", class = "Intercept", nlpar = "a")
  expect_identical(as.data.frame(bi)[, c("class", "coef", "nlpar")],
                   data.frame(class = "b", coef = "Intercept", nlpar = "a"))
  f <- bf(yn ~ a * exp(k * x), a ~ 1, k ~ 1, nl = TRUE)
  expect_error(
    frm(f, data = d, family = gaussian(),
        prior = set_prior("normal(7, 1)", nlpar = "a", coef = "Intercept") +
          bi + set_prior("normal(0, 1)", nlpar = "k")),
    "Duplicated prior specifications", class = "frmtmb_error")
})

test_that("no spelling left reaches a parameter a class reaches", {
  # the three spellings that once did: class theta with a coef, class
  # theta with none, and a named list of prior objects keyed by internal
  # name (which also dropped the LKJ on the rest of the block)
  d <- ros_data()
  for (f in list(bf(y ~ x + z + (1 + x + z | g)),
                 bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE))) {
    des <- frmtmb:::prior_design(f, d, NULL, list())
    cc <- frmtmb:::cov_class_coefs(des$spec, des$frame)
    expect_gt(nrow(cc), 0L)
    for (nm in cc$name) {
      expect_error(ros_reached(f, d, set_prior("normal(0, 1)",
                                               class = "theta",
                                               coef = nm)),
                   "has a class of its own", class = "frmtmb_error")
      pl <- list(prior_normal(0, 1))
      names(pl) <- nm
      expect_error(frmtmb:::resolve_prior_input(
        frmtmb:::prior_design(f, d, NULL, list()), pl),
        "no longer takes it", class = "frmtmb_error")
    }
    expect_error(ros_reached(f, d, set_prior("normal(0, 1)",
                                             class = "theta")),
                 "names no parameter", class = "frmtmb_error")
  }
})

test_that("\"\" names a slot and puts nothing on it in frm()", {
  # brms's (flat) on one slot: frm() reads it as no prior, so a fit
  # given only such slots is the maximum-likelihood fit
  d <- ros_data()
  f <- bf(y ~ x + (1 + x | g))
  ml <- frm(f, data = d, family = gaussian())
  for (pr in list(set_prior("", class = "cor"),
                  set_prior("", class = "sd", group = "g"),
                  set_prior("", class = "Intercept"))) {
    fl <- frm(f, data = d, family = gaussian(), prior = pr)
    expect_null(fl$prior)
    expect_identical(fl$opt$par, ml$opt$par)
  }
  expect_length(ros_reached(f, d, set_prior("", class = "cor")), 0L)
  # and the class still has to exist
  expect_error(validate_prior(set_prior("", class = "rescor"), f,
                              data = d),
               "No residual correlation matches", class = "frmtmb_error")
})

test_that("a cor refusal with resp says that resp is frmtmb's own", {
  d <- ros_data()
  f <- bf(y ~ x + (1 + x | g)) + bf(y2 ~ x + (1 + x | g)) +
    set_rescor(FALSE)
  m <- ros_msg(validate_prior(set_prior("normal(0, 1)", class = "theta",
                                        coef = "theta_3"), f, data = d))
  expect_match(m, "class = \"cor\", group = \"g\", resp = \"y\"",
               fixed = TRUE)
  expect_match(m, "frmtmb's own spelling: brms refuses it", fixed = TRUE)
})

test_that("Intercept with nlpar is named as written when it misses", {
  d <- ros_data()
  d$yn <- 2 * exp(0.3 * d$x) + stats::rnorm(nrow(d), 0, 0.3)
  m <- ros_msg(validate_prior(
    set_prior("normal(0, 1)", class = "Intercept", nlpar = "a"),
    bf(yn ~ a * exp(k * x), a ~ 0 + x, k ~ 1, nl = TRUE), data = d))
  expect_match(m, "Prior target not found (class=Intercept, nlpar=a)",
               fixed = TRUE)
})

test_that("c() and a list refuse what is not a specification", {
  sp <- set_prior("normal(0, 1)", class = "b")
  expect_error(c(sp, list(theta_4 = prior_normal(0, 1))),
               "prior = list(theta_4 = ) put densities on internal",
               fixed = TRUE)
  expect_error(c(sp, prior_normal(0, 1)),
               "a density object alone names no parameter", fixed = TRUE)
  d <- ros_data()
  expect_error(frm(bf(y ~ x + (1 | g)), data = d,
                   prior = list(sp, theta_4 = prior_normal(0, 1))),
               "prior = list(theta_4 = ) put densities on internal",
               fixed = TRUE)
})

test_that("frm()'s table carries its route into a priorlist", {
  # frmtmb.sample refuses it by this mark; frm() reads the table, and
  # the prior it keeps carries no mark
  d <- ros_data()
  f <- y ~ x + (1 + x | g)
  p <- default_prior(f, data = d)
  p$prior[p$class == "b" & p$coef == "x"] <- "normal(0, 1)"
  sp <- set_prior("normal(0, 1)", class = "sigma")
  for (pl in list(frmtmb:::as_priorlist(p), as.brmsprior(p),
                  sp + as.brmsprior(p), as.brmsprior(p) + sp,
                  c(sp, as.brmsprior(p)))) {
    expect_identical(attr(pl, "route"), "fit")
  }
  fit <- frm(f, data = d, family = gaussian(), prior = p)
  expect_null(attr(fit$prior, "route"))
  # and the fit is the one the edited row alone gives
  fx <- frm(f, data = d, family = gaussian(),
            prior = set_prior("normal(0, 1)", class = "b", coef = "x"))
  expect_identical(fit$opt$par, fx$opt$par)
})

test_that("a list of specifications points to + and c()", {
  # the 0.77.0 consolidation: list(sp, sp2) got the generic "not a list
  # of length 2"; brms combines specifications with + or c()
  d <- ros_data()
  sp <- set_prior("normal(0, 1)", class = "b")
  sp2 <- set_prior("normal(0, 2)", class = "Intercept")
  expect_error(frm(bf(y ~ x), data = d, prior = list(sp, sp2)),
               "Combine them with + or c(), as in brms: sp1 + sp2",
               fixed = TRUE)
  # the guard absent: both spellings fit, and give the same fit
  f1 <- frm(bf(y ~ x), data = d, prior = sp + sp2)
  f2 <- frm(bf(y ~ x), data = d, prior = c(sp, sp2))
  expect_identical(f1$opt$par, f2$opt$par)
})

test_that("c() reads a prior table as brms's c() does, keeping its route", {
  # the 0.77.0 release review, finding 3: c(set_prior(), default_prior())
  # was refused as "a plain data frame"; 0.76.0 combined them, and brms
  # 2.23.0's c() combines a brmsprior and a prior table in either order
  # (dev/rel077-fix2-brms-c.R)
  d <- ros_data()
  f <- y ~ x + (1 + x | g)
  pf <- default_prior(f, data = d)
  sp <- set_prior("normal(0, 1)", class = "b", coef = "x")
  for (pl in list(c(sp, pf), c(pf, sp))) {
    expect_s3_class(pl, "frmtmb_priorlist")
    expect_identical(attr(pl, "route"), "fit")
    expect_length(pl, length(as.brmsprior(pf)) + 1L)
  }
  # the rows are the specification and the table read as as.brmsprior()
  # reads it
  expect_equal(as.data.frame(c(sp, pf))[-1L, ],
               as.data.frame(as.brmsprior(pf)), ignore_attr = TRUE)
  # the guard absent: what is not a specification or a table is refused
  expect_error(c(sp, prior_normal(0, 1)),
               "a density object alone names no parameter", fixed = TRUE)
})
