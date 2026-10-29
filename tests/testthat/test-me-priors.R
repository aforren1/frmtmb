# brms's prior classes of the me() hyperparameters: meanme, sdme and
# corme (dev/meme-findings.md). brms 2.23.0's rows, bounds and refusals
# are in dev/meme-out/brms.txt (script dev/meme-brms.R). Through 0.77.0
# frmtmb listed no row for the three classes and refused them by name.

mp_data <- function(seed = 51, n = 60) {
  set.seed(seed)
  tx <- rnorm(n, 1, 0.8)
  tz <- 0.5 * tx + rnorm(n, 0, 0.7)
  d <- data.frame(x = tx + rnorm(n, 0, 0.3), sx = 0.3,
                  z = tz + rnorm(n, 0, 0.3), sz = runif(n, 0.2, 0.4),
                  g = factor(rep(seq_len(n / 5), each = 5)),
                  h = factor(rep(1:5, n / 5)))
  ng <- nlevels(d$g)
  d$xg <- rnorm(ng)[d$g]
  d$zg <- 0.6 * d$xg + rnorm(ng, 0, 0.6)[d$g]
  d$sg <- 0.25
  d$y <- 2 + 0.7 * tx - 0.4 * tz + 0.3 * d$xg + rnorm(n, 0, 0.5)
  d
}

# the me() rows of a prior table as class|coef|group|lb
mp_rows <- function(tab) {
  tab <- as.data.frame(tab)
  tab <- tab[tab$class %in% c("meanme", "sdme", "corme"), , drop = FALSE]
  lb <- suppressWarnings(as.numeric(tab$lb))
  paste(tab$class, tab$coef, tab$group, ifelse(is.na(lb), "", lb),
        sep = "|")
}

mp_forms <- list(
  one = y ~ me(x, sx),
  two = y ~ me(x, sx) + me(z, sz),
  nocor = bf(y ~ me(x, sx) + me(z, sz)) + set_mecor(FALSE),
  grouped = y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g),
  mixed = y ~ me(x, sx) + me(z, sz) + me(xg, sg, gr = g))

test_that("default_prior() lists brms's meanme, sdme and corme rows", {
  d <- mp_data()
  # brms 2.23.0's rows for these shapes, in its order
  # (dev/meme-out/brms.txt, cases A, B, C, E and G)
  want <- list(
    one = c("meanme|||", "meanme|mex||", "sdme|||0", "sdme|mex||"),
    two = c("corme|||", "meanme|||", "meanme|mex||", "meanme|mez||",
            "sdme|||0", "sdme|mex||", "sdme|mez||"),
    nocor = c("meanme|||", "meanme|mex||", "meanme|mez||", "sdme|||0",
              "sdme|mex||", "sdme|mez||"),
    grouped = c("corme|||", "corme||g|", "meanme|||", "meanme|mexg||",
                "meanme|mezg||", "sdme|||0", "sdme|mexg||", "sdme|mezg||"),
    mixed = c("corme|||", "meanme|||", "meanme|mex||", "meanme|mexg||",
              "meanme|mez||", "sdme|||0", "sdme|mex||", "sdme|mexg||",
              "sdme|mez||"))
  for (nm in names(mp_forms)) {
    tab <- default_prior(mp_forms[[nm]], data = d)
    expect_identical(mp_rows(tab), want[[nm]], label = nm)
    # frm() applies no density on these slots
    me_tab <- as.data.frame(tab)
    me_tab <- me_tab[me_tab$class %in% c("meanme", "sdme", "corme"), ]
    expect_true(all(me_tab$prior == "(flat)"), label = nm)
  }
})

test_that("the me() rows are brms's own, read from brms", {
  skip_unless_brms()
  d <- mp_data()
  brms_forms <- list(
    one = brms::bf(y ~ me(x, sx)),
    two = brms::bf(y ~ me(x, sx) + me(z, sz)),
    nocor = brms::bf(y ~ me(x, sx) + me(z, sz)) + brms::set_mecor(FALSE),
    grouped = brms::bf(y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g)),
    mixed = brms::bf(y ~ me(x, sx) + me(z, sz) + me(xg, sg, gr = g)))
  for (nm in names(brms_forms)) {
    bt <- with_brms_me(function() {
      brms::default_prior(brms_forms[[nm]], data = d)
    })
    ft <- default_prior(mp_forms[[nm]], data = d)
    expect_identical(mp_rows(ft), mp_rows(bt), label = nm)
  }
})

test_that("a model without me() has no me() rows and refuses the classes", {
  d <- mp_data()
  tab <- as.data.frame(default_prior(y ~ x + (1 | h), data = d))
  expect_false(any(tab$class %in% c("meanme", "sdme", "corme")))
  for (p in list(set_prior("normal(0, 1)", class = "meanme"),
                 set_prior("exponential(1)", class = "sdme"),
                 set_prior("lkj(2)", class = "corme"))) {
    expect_error(frm(y ~ x + (1 | h), data = d, prior = p),
                 "this model has no me[(][)] term. brms refuses the row too")
  }
})

test_that("set_prior() refuses what brms's me() rows cannot carry", {
  expect_error(set_prior("normal(0, 1)", class = "meanme", group = "g"),
               "class = \"meanme\" has no group")
  expect_error(set_prior("exponential(1)", class = "sdme", group = "g"),
               "class = \"sdme\" has no group")
  expect_error(set_prior("lkj(2)", class = "corme", coef = "mex"),
               "names nothing it can narrow to")
  # brms: "The only supported prior for correlation matrices is the
  # 'lkj' prior"
  expect_error(set_prior("normal(0, 1)", class = "corme"),
               "class = \"corme\" takes an lkj[(][)] prior")
  expect_error(set_prior("lkj(2)", class = "corme", lb = 0),
               "class = \"corme\" takes no lb/ub")
  for (a in list(list(resp = "y"), list(dpar = "sigma"),
                 list(nlpar = "a"))) {
    expect_error(do.call(set_prior, c(list("normal(0, 1)",
                                           class = "meanme"), a)),
                 paste0("class = \"meanme\" takes no ", names(a)))
  }
  # lkj() belongs to the correlation classes only
  expect_error(set_prior("lkj(2)", class = "sdme"), "belongs to class")
})

test_that("a prior names what the model has, or is refused naming it", {
  d <- mp_data()
  expect_error(frm(mp_forms$two, data = d,
                   prior = set_prior("normal(0, 1)", class = "meanme",
                                     coef = "mew")),
               "this model's are \"mex\", \"mez\"")
  expect_error(frm(mp_forms$nocor, data = d,
                   prior = set_prior("lkj(2)", class = "corme")),
               "set_mecor[(]FALSE[)] fixes it at zero")
  expect_error(frm(mp_forms$one, data = d,
                   prior = set_prior("lkj(2)", class = "corme")),
               "and this model has none")
  expect_error(frm(mp_forms$grouped, data = d,
                   prior = set_prior("lkj(2)", class = "corme",
                                     group = "h")),
               "grouped by \"g\"")
  expect_error(frm(mp_forms$two, data = d,
                   prior = set_prior("lkj(2)", class = "corme",
                                     group = "g")),
               "have no gr =")
})

test_that("each class resolves onto its own parameter, one slot each", {
  d <- mp_data()
  design <- frmtmb:::prior_design(mp_forms$mixed, d, NULL, list())
  pr <-set_prior("normal(0, 2)", class = "meanme") +
    set_prior("normal(1, 0.5)", class = "meanme", coef = "mez") +
    set_prior("exponential(1)", class = "sdme") +
    set_prior("student_t(3, 0, 1)", class = "sdme", coef = "mex") +
    set_prior("lkj(2)", class = "corme")
  ri <- frmtmb:::resolve_priorlist(design, pr)
  key <- vapply(ri$entries, function(e) {
    paste0(e$comp, ".", paste(e$idx, collapse = ","))
  }, "")
  # three meanme, three logsdme and the one correlated pair
  expect_setequal(key, c("meanme.1", "meanme.2", "meanme.3", "logsdme.1",
                         "logsdme.2", "logsdme.3", "thetame.1"))
  at <- function(k) ri$entries[[match(k, key)]]
  expect_identical(at("meanme.2")$dist$scale, 0.5)
  expect_identical(at("meanme.1")$dist$scale, 2)
  expect_identical(at("meanme.1")$scale, "internal")
  expect_identical(at("logsdme.1")$dist$kind, "t")
  expect_identical(at("logsdme.2")$dist$kind, "exponential")
  expect_identical(at("logsdme.1")$scale, "sd")
  expect_identical(at("thetame.1")$dist$kind, "lkj")
  expect_identical(at("thetame.1")$dist$map, "chol")
  # the coef row wins over the class row whatever the written order
  pr2 <- set_prior("normal(1, 0.5)", class = "meanme", coef = "mez") +
    set_prior("normal(0, 2)", class = "meanme")
  ri2 <- frmtmb:::resolve_priorlist(design, pr2)
  k2 <- vapply(ri2$entries, function(e) paste0(e$comp, ".", e$idx), "")
  expect_identical(ri2$entries[[match("meanme.2", k2)]]$dist$scale, 0.5)
  # a bound is on the parameter itself: sdme's on its log
  pr3 <- set_prior("", class = "sdme", lb = 0.5, ub = 4) +
    set_prior("", class = "meanme", coef = "mex", lb = -1)
  ri3 <- frmtmb:::resolve_priorlist(design, pr3)
  expect_identical(unname(ri3$lower[c("logsdme_mex", "logsdme_mez",
                                      "logsdme_mexg")]), rep(log(0.5), 3))
  expect_identical(unname(ri3$upper[["logsdme_mex"]]), log(4))
  expect_identical(unname(ri3$lower[["meanme_mex"]]), -1)
  expect_length(ri3$entries, 0L)
})

test_that("a group's corme row replaces the class row on that group", {
  d <- mp_data(seed = 52)
  d$wg <- rnorm(nlevels(d$g))[d$g]
  form <- y ~ me(x, sx) + me(z, sz) + me(xg, sg, gr = g) + me(wg, sg, gr = g)
  design <- frmtmb:::prior_design(form, d, NULL, list())
  me <- design$frame$me
  expect_identical(vapply(me$groups, `[[`, "", "gr"), c("", "g"))
  pr <- set_prior("lkj(2)", class = "corme") +
    set_prior("lkj(4)", class = "corme", group = "g")
  ri <- frmtmb:::resolve_priorlist(design, pr)
  eta <- vapply(ri$entries, function(e) e$dist$eta, 0)
  idx <- lapply(ri$entries, `[[`, "idx")
  expect_identical(eta[[match(list(me$groups[[1L]]$th_idx), idx)]], 2)
  expect_identical(eta[[match(list(me$groups[[2L]]$th_idx), idx)]], 4)
  # written in the other order, the same
  ri2 <- frmtmb:::resolve_priorlist(design, pr[2:1])
  expect_setequal(vapply(ri2$entries, function(e) e$dist$eta, 0), c(2, 4))
  # a table shows it
  vp <- as.data.frame(validate_prior(pr, form, data = d))
  vp <- vp[vp$class == "corme", ]
  expect_identical(vp$prior, c("lkj(2)", "lkj(4)"))
  expect_identical(vp$group, c("", "g"))
})

# the log density frmtmb adds for these priors at one point, written
# from the densities themselves: brms's lpdf of each statement, the log
# map's Jacobian for sdme, and for a pair the LKJ(eta) density of the
# correlation, (1 - rho^2)^(eta - 1) / (2^(2 eta - 1) B(eta, eta)), with
# the Jacobian of rho = t / sqrt(1 + t^2)
mp_ref_logprior <- function(meanme, logsdme, t, eta) {
  s <- exp(logsdme)
  rho <- t / sqrt(1 + t^2)
  stats::dnorm(meanme[[1]], 0, 2, log = TRUE) +
    stats::dnorm(meanme[[2]], 1, 0.5, log = TRUE) +
    stats::dt(s[[1]], 3, log = TRUE) + logsdme[[1]] +
    stats::dexp(s[[2]], 1, log = TRUE) + logsdme[[2]] +
    (eta - 1) * log(1 - rho^2) - (2 * eta - 1) * log(2) -
    lbeta(eta, eta) - 1.5 * log(1 + t^2)
}

test_that("frm() adds brms's densities to the objective, and moves the mode", {
  d <- mp_data()
  pr <- set_prior("normal(0, 2)", class = "meanme") +
    set_prior("normal(1, 0.5)", class = "meanme", coef = "mez") +
    set_prior("exponential(1)", class = "sdme") +
    set_prior("student_t(3, 0, 1)", class = "sdme", coef = "mex") +
    set_prior("lkj(2)", class = "corme")
  f0 <- frm(mp_forms$two, data = d)
  fp <- frm(mp_forms$two, data = d, prior = pr)
  expect_identical(names(fp$obj$env$last.par.best),
                   names(f0$obj$env$last.par.best))
  for (p in list(fp$obj$env$last.par.best, f0$obj$env$last.par.best)) {
    nm <- names(p)
    diff <- f0$obj$env$f(p) - fp$obj$env$f(p)
    ref <- mp_ref_logprior(p[nm == "meanme"], p[nm == "logsdme"],
                           p[nm == "thetame"], 2)
    expect_lt(abs(diff - ref), 1e-10 * max(1, abs(ref)))
  }
  # a tight prior pulls its parameter to the prior's location
  ft <- frm(mp_forms$two, data = d,
            prior = set_prior("normal(5, 0.001)", class = "meanme",
                              coef = "mex"))
  m0 <- f0$estimates$meanme[["meanme_mex"]]
  mt <- ft$estimates$meanme[["meanme_mex"]]
  expect_lt(abs(mt - 5), abs(m0 - 5) / 100)
  # prior_summary() and the method line name the penalty
  ps <- as.data.frame(prior_summary(fp))
  expect_setequal(paste(ps$class, ps$coef),
                  c("meanme ", "meanme mez", "sdme ", "sdme mex", "corme "))
})

test_that("validate_prior() fills the me() rows as it fills class b's", {
  # a coefficient row under a user class row reads "(vectorized)", as
  # brms PRINTS its table and as frmtmb's class "b" rows read; brms's
  # data frame holds "" and "default" there
  d <- mp_data()
  pr <- set_prior("normal(0, 2)", class = "meanme") +
    set_prior("normal(1, 0.5)", class = "meanme", coef = "mez") +
    set_prior("exponential(1)", class = "sdme")
  vp <- as.data.frame(validate_prior(pr, mp_forms$two, data = d))
  vp <- vp[vp$class %in% c("meanme", "sdme"), ]
  got <- stats::setNames(paste(vp$prior, vp$source),
                         paste(vp$class, vp$coef))
  expect_identical(got[["meanme "]], "normal(0, 2) user")
  expect_identical(got[["meanme mez"]], "normal(1, 0.5) user")
  expect_identical(got[["meanme mex"]], "normal(0, 2) (vectorized)")
  expect_identical(got[["sdme "]], "exponential(1) user")
  expect_identical(got[["sdme mex"]], "exponential(1) (vectorized)")
  # a table passed back is the same prior
  f1 <- frm(mp_forms$two, data = d, prior = pr)
  f2 <- frm(mp_forms$two, data = d,
            prior = validate_prior(pr, mp_forms$two, data = d))
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f2)))
})

test_that("default_prior() answers me(x) written without sdx; frm() refuses", {
  d <- mp_data()
  # brms 2.23.0's default_prior() lists these rows and only standata()
  # refuses (dev/meme-out/fix-brms-grp.txt); through 0.77.0 frmtmb
  # refused at the parse, so no table answered
  expect_identical(mp_rows(default_prior(y ~ me(x), data = d)),
                   c("meanme|||", "meanme|mex||", "sdme|||0", "sdme|mex||"))
  vp <- as.data.frame(validate_prior(
    set_prior("normal(0, 1)", class = "meanme"), y ~ me(x), data = d))
  expect_identical(vp$prior[vp$class == "meanme" & !nzchar(vp$coef)],
                   "normal(0, 1)")
  expect_error(frm(y ~ me(x), data = d),
               "Argument 'sdx' is missing in function 'me': me[(]x[)]")
})
