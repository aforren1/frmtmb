## brms's link arguments for wiener(): which links each dpar takes, the
## argument order, and what a non-default link must and must not change.

test_that("wiener() takes brms's link arguments, in brms's order", {
  # brms's wiener(link, link_bs, link_ndt, link_bias), so that a
  # positional call written for brms means the same thing here
  expect_equal(names(formals(wiener))[1:4],
               c("link", "link_bs", "link_ndt", "link_bias"))
  skip_if_not_installed("brms")
  expect_equal(names(formals(wiener))[1:4], names(formals(brms::wiener)))
})

test_that("the mean link is refused outside brms's set", {
  # Before, wiener(link = "logit") FITTED: a drift rate squeezed into
  # (0, 1), which brms refuses.
  for (lk in c("logit", "probit", "inverse", "sqrt", "cloglog")) {
    expect_error(wiener(link = lk),
                 "is not a supported link for family 'wiener'",
                 class = "frmtmb_eam_error")
  }
  for (lk in ddm_links_brms) {
    expect_equal(wiener(link = lk)$links$mu$name, lk)
  }
  skip_if_not_installed("brms")
  expect_setequal(ddm_links_brms,
                  get(".family_wiener", asNamespace("brms"))()$links)
})

test_that("an unquoted link name is read as its name, as brms reads it", {
  expect_equal(wiener(log)$links$mu$name, "log")
  expect_equal(wiener(softplus)$links$mu$name, "softplus")
  # a variable holding a name is still evaluated
  lk <- "squareplus"
  expect_equal(wiener(lk)$links$mu$name, "squareplus")
})

test_that("link_bs, link_ndt and link_bias are refused outside brms's sets", {
  expect_error(wiener(link_bs = "logit"),
               "is not a supported link for parameter 'bs'",
               class = "frmtmb_eam_error")
  expect_error(wiener(link_bias = "log"),
               "is not a supported link for parameter 'bias'",
               class = "frmtmb_eam_error")
  expect_error(wiener(link_ndt = "logit"),
               "is not a supported link for parameter 'ndt'",
               class = "frmtmb_eam_error")
  expect_error(wiener(link_ndt = c("log", "identity")),
               "takes a single link name", class = "frmtmb_eam_error")
})

test_that("a link on a scalar dpar moves neither maximum nor estimate", {
  # A bijective link on a dpar with no formula reparameterizes the same
  # model, so the maximized log-likelihood and the natural-scale
  # estimate must come back as the default link's, to the optimizer's
  # own resolution, which the default fit's standard errors measure.
  set.seed(20261013)
  d <- ddm_simulate(300, mu = 0.8, bs = 1.3, ndt = 0.3, bias = 0.45)
  d$x <- stats::rnorm(300)
  f0 <- frm(bf(rt | dec(upper) ~ x), family = wiener(), data = d)
  ll0 <- as.numeric(logLik(f0))
  sp0 <- summary(f0)$spec_pars
  arms <- list(
    bs = list(fam = wiener(link_bs = "softplus"), dpar = "bs",
              link = "softplus"),
    ndt = list(fam = wiener(link_ndt = "log"), dpar = "ndt", link = "log"),
    ndt_id = list(fam = wiener(link_ndt = "identity"), dpar = "ndt",
                  link = "identity"),
    bias = list(fam = wiener(link_bias = "identity"), dpar = "bias",
                link = "identity"))
  for (nm in names(arms)) {
    a <- arms[[nm]]
    f <- frm(bf(rt | dec(upper) ~ x), family = a$fam, data = d)
    expect_equal(f$opt$convergence, 0L, info = nm)
    expect_equal(family(f)$links[[a$dpar]]$name, a$link, info = nm)
    expect_lt(abs(as.numeric(logLik(f)) - ll0), 1e-8 * abs(ll0))
    sp <- summary(f)$spec_pars
    expect_lt(max(abs(sp[, "Estimate"] - sp0[rownames(sp), "Estimate"]) /
                    sp0[rownames(sp), "Est.Error"]), 1e-3)
  }
})

test_that("link_ndt = 'log' is brms's parameterization of ndt ~ formula", {
  set.seed(20261014)
  d <- ddm_simulate(300, mu = 0.8, bs = 1.3, ndt = 0.3, bias = 0.5)
  d$x <- stats::rnorm(300)
  f <- frm(bf(rt | dec(upper) ~ 1, ndt ~ x), family = wiener(link_ndt = "log"),
           data = d)
  expect_equal(f$opt$convergence, 0L)
  # the coefficients are on the log scale, as brms's are
  eta <- frmtmb::frm_linpred(f, dpar = "ndt", type = "link")
  t0 <- frmtmb::frm_linpred(f, dpar = "ndt", type = "response")
  expect_equal(as.numeric(t0), exp(as.numeric(eta)))
  expect_equal(as.numeric(ndt_time(f)), as.numeric(t0))
  # each row's own density keeps its non-decision time below its time
  expect_true(all(as.numeric(t0) < d$rt))
  skip_if_not_installed("brms")
  expect_equal(brms::wiener()$link_ndt, "log")
  sc <- as.character(brms::stancode(
    brms::bf(rt | dec(upper) ~ 1, ndt ~ x), family = brms::wiener(),
    data = d))
  expect_match(sc, "ndt = exp(ndt)", fixed = TRUE)
})

test_that("link_ndt is refused with whatever needs the bound", {
  set.seed(20261015)
  d <- ddm_simulate(80, mu = 0.8, bs = 1.3, ndt = 0.3)
  d$s <- factor(rep(1:4, each = 20))
  expect_error(wiener(link_ndt = "log", max_ndt = 0.2),
               "`max_ndt` is the bound of the default",
               class = "frmtmb_eam_error")
  expect_error(frm(bf(rt | dec(upper) + ndt_group(s) ~ 1),
                   family = wiener(link_ndt = "log"), data = d),
               "ndt_group\\(\\) sets the bound", class = "frmtmb_eam_error")
  # every row pinned at the floored decision time at once is a spike
  # the likelihood has no upper bound on (dev/famuni-wiener-findings.md)
  for (a in list(list(variability = "sv"), list(contaminant = TRUE),
                 list(allow_unreachable = TRUE))) {
    expect_error(do.call(wiener, c(list(link_ndt = "log"), a)),
                 "brms's parameterization of the plain model",
                 class = "frmtmb_eam_error")
  }
  # and the default link keeps every one of them
  expect_no_error(wiener(max_ndt = 0.2))
  expect_no_error(wiener(variability = "sv"))
  expect_no_error(wiener(max_ndt = 0.4, allow_unreachable = TRUE))
})

test_that("an identity link returns NaN outside the support, not a number", {
  # ddm_pos() maps a negative series sum to its magnitude, so the plain
  # density at a bias of 1.2 is a finite number that is not a density.
  # Only the default links can never reach there.
  set.seed(20261016)
  d <- ddm_simulate(60, mu = 0.8, bs = 1.3, ndt = 0.3)
  fid <- frm(bf(rt | dec(upper) ~ 1),
             family = wiener(link_bs = "identity", link_ndt = "identity",
                             link_bias = "identity"), data = d)
  f0 <- frm(bf(rt | dec(upper) ~ 1), family = wiener(), data = d)
  y <- d$rt[1:3]
  av <- list(dec = d$upper[1:3])
  at <- function(bs = 1.3, bias = 0.45, ndt = 0.2) {
    list(mu = 0.8, bs = bs, ndt = ndt, bias = bias)
  }
  inside <- family(fid)$lpdf(y, at(), av)
  expect_true(all(is.finite(inside)))
  expect_identical(inside, family(f0)$lpdf(y, at(), av))
  expect_true(all(is.finite(family(f0)$lpdf(y, at(bias = 1.2), av))))
  expect_true(all(is.nan(family(fid)$lpdf(y, at(bias = 1.2), av))))
  expect_true(all(is.nan(family(fid)$lpdf(y, at(ndt = -0.05), av))))
  allow_warnings(
    expect_true(all(is.nan(family(fid)$lpdf(y, at(bs = -1), av)))),
    "NaNs produced")
})

test_that("the fitted objective has no value outside the support", {
  # The test above reads the density on doubles. A fit reads it on the
  # RTMB tape, which folds `0 * x` and `0 / x` to a constant 0 when it
  # records them: the first guard, 0 / (x + |x|), was NaN on doubles and
  # 0 on the tape, so at a bias of -0.5 the objective was finite and its
  # log-likelihood (+10.57) above the maximum (dev/reviews/
  # 2026-10-13-famuni.md, fix 1). This reads the objective itself.
  set.seed(20261016)
  d <- ddm_simulate(60, mu = 0.8, bs = 1.3, ndt = 0.3)
  fid <- frm(bf(rt | dec(upper) ~ 1),
             family = wiener(link_bs = "identity", link_ndt = "identity",
                             link_bias = "identity"), data = d)
  par <- fid$opt$par
  expect_true(is.finite(fid$obj$fn(par)))
  # betad holds bs, ndt and bias, in that order, each on its own scale
  ib <- which(names(par) == "betad")
  expect_length(ib, 3L)
  for (k in seq_along(ib)) {
    for (v in c(-0.5, if (k == 3L) 1.2)) {
      p2 <- par
      p2[ib[k]] <- v
      expect_true(is.nan(fid$obj$fn(p2)),
                  label = paste0("objective at betad[", k, "] = ", v))
    }
  }
})
