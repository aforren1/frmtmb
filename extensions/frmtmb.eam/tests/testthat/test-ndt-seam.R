## The SECOND thing this package promises another extension package: the
## bound a non-decision time is measured against. R/ndt-seam.R says what
## is and is not promised; this file is that promise, asserted.
##
## Item 1.0b of dev/extension-gaps-plan.md asked for the seam because
## frmtmb.learn::rlddm() had written its own copy of the scaled logit
## and so inherited the one-global-bound defect item 1.0a removed here.
## What the seam has to be worth is that the copy goes away and the
## behavior does not.

test_that("ndt_bound derives one bound, or one per group", {
  rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
  g <- c("a", "a", "a", "b", "b")
  b0 <- ndt_bound(rt)
  expect_s3_class(b0, "frmtmb_eam_ndt_bound")
  expect_identical(b0$ub, min(rt))
  expect_null(b0$floors)
  expect_false(isTRUE(b0$pending))

  bg <- ndt_bound(rt, list(ndt_group = ndt_bound_key(g)))
  expect_identical(sort(unname(bg$floors)), c(0.31, 0.61))
  expect_identical(sort(unname(bg$sizes)), c(2L, 3L))
  # the table is keyed on the LABEL's code, so it is the same table
  # whatever else was in the column the grouping was read from
  expect_identical(bg$floors,
                   ndt_bound(rt, list(ndt_group = ndt_bound_key(
                     factor(g, levels = c("b", "a")))))$floors)

  # max_ndt is the stated bound, and above the fastest response it is
  # refused rather than accepted
  expect_identical(ndt_bound(rt, max_ndt = 0.2)$ub, 0.2)
  expect_error(ndt_bound(rt, max_ndt = 5, what = "my_family"),
               "my_family: max_ndt")
  expect_error(ndt_bound(rt, list(ndt_group = ndt_bound_key(g)),
                         max_ndt = 0.2),
               "both set the non-decision time's upper bound")
})

test_that("ndt_bound_key is a function of the label and of nothing else", {
  g <- c("s3", "s1", "s2", "s1")
  a <- ndt_bound_key(g)
  expect_identical(a, ndt_bound_key(factor(g)))
  expect_identical(a, ndt_bound_key(factor(g, levels = c("s3", "s2", "s1"))))
  expect_identical(a, ndt_bound_key(droplevels(factor(g)[1:3]))[c(1:3, 2)])
  expect_identical(ndt_bound_key(c(TRUE, FALSE)),
                   ndt_bound_key(c("TRUE", "FALSE")))
  expect_error(ndt_bound_key(c("a", NA)), "every row needs one")
})

test_that("ndt_bound_pending is the state a family has before frm()", {
  p <- ndt_bound_pending(what = "my_family")
  expect_true(isTRUE(p$pending))
  expect_null(ndt_bound_of(list(ndt_bound = p)))
  s <- ndt_bound_pending(0.4, what = "my_family")
  expect_false(isTRUE(s$pending))
  expect_identical(s$ub, 0.4)
  expect_identical(ndt_bound_of(list(ndt_bound = s))$ub, 0.4)
  expect_error(ndt_bound_pending(c(1, 2), what = "my_family"),
               "one positive number")
  expect_error(ndt_bound_pending(-1, what = "my_family"),
               "one positive number")
})

# A family of another package's shape: a plain frmtmb family with an
# `ndt` parameter and no relation to this package's four. It is the
# whole point of the seam that it works on one of these, so the test
# builds one rather than borrowing wiener().
seam_family <- function() {
  frmtmb::frmtmb_family(
    "seamtest",
    dpars = c("mu", "ndt"),
    links = list(mu = "identity", ndt = "log"),
    primary_dpars = "mu",
    type = "continuous",
    # the seam's worked example: one call where it would have read
    # dpars[["ndt"]], and nothing else about the density changes
    lpdf = function(y, dpars, aterms) {
      stats::dnorm(y - ndt_apply(dpars, aterms, "seamtest"),
                   dpars[["mu"]], 1, log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) 0,
                      ndt = function(y, aterms) 0.5 * min(y)))
}

test_that("ndt_bound_attach puts a scalar bound in the link and nothing else", {
  rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
  fam <- seam_family()
  init0 <- fam$init_dpars$ndt
  out <- ndt_bound_attach(fam, ndt_bound(rt, what = "seamtest"))
  expect_identical(out$links$ndt$name, "scaled_logit")
  # the bound is IN the link, so `ndt` is a time on the response scale
  expect_equal(out$links$ndt$linkinv(0), min(rt) / 2)
  expect_equal(out$links$ndt$linkfun(min(rt) / 2), 0)
  # and no per-row data is added, because there is no per-row bound
  expect_null(out$aterm_data)
  # the family's own starting value is left alone
  expect_identical(out$init_dpars$ndt, init0)
  expect_identical(ndt_bound_of(out)$ub, min(rt))
})

test_that("a per-group bound becomes a logit plus per-row data", {
  rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
  g <- c("a", "a", "a", "b", "b")
  at <- list(ndt_group = ndt_bound_key(g))
  out <- ndt_bound_attach(seam_family(), ndt_bound(rt, at, what = "s"))
  expect_identical(out$links$ndt, "logit")
  expect_equal(out$init_dpars$ndt(rt, at), 0.5)
  av <- out$aterm_data(rt, at)
  expect_named(av, "ndt_floor")
  expect_equal(av$ndt_floor, c(0.31, 0.31, 0.31, 0.61, 0.61))
  # the floor is read out of the CAPTURED table, not recomputed from
  # whatever response this is handed: a leave-one-out refit that dropped
  # a group's fastest trial must score against the bound the fit used
  expect_equal(out$aterm_data(rt + 10, at)$ndt_floor,
               c(0.31, 0.31, 0.31, 0.61, 0.61))
})

test_that("attaching twice replaces the bound instead of stacking on it", {
  # frm() runs family_finalize() again on the influence(),
  # frm_simulate() and prior-predictive paths, and the failure this
  # guards is a SECOND ndt_floor entry in the addition-term values.
  rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
  at <- list(ndt_group = ndt_bound_key(c("a", "a", "a", "b", "b")))
  one <- ndt_bound_attach(seam_family(), ndt_bound(rt, at, what = "s"))
  two <- ndt_bound_attach(one, ndt_bound(rt, at, what = "s"))
  expect_length(two$aterm_data(rt, at), 1L)
  expect_equal(two$aterm_data(rt, at)$ndt_floor,
               one$aterm_data(rt, at)$ndt_floor)
  # and a scalar bound attached over a per-group one leaves no per-row
  # data behind at all
  back <- ndt_bound_attach(one, ndt_bound(rt, what = "s"))
  expect_null(back$aterm_data)
  expect_identical(back$links$ndt$name, "scaled_logit")
})

test_that("a family that composes its own aterm_data keeps it", {
  rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
  at <- list(ndt_group = ndt_bound_key(c("a", "a", "a", "b", "b")))
  fam <- seam_family()
  fam$aterm_data <- function(y, aterms) list(own = seq_along(y))
  out <- ndt_bound_attach(fam, ndt_bound(rt, at, what = "s"))
  expect_named(out$aterm_data(rt, at), c("own", "ndt_floor"))
})

test_that("ndt_bound_attach refuses what it cannot honestly do", {
  rt <- c(0.31, 0.42, 0.55)
  b <- ndt_bound(rt, what = "s")
  expect_error(ndt_bound_attach(list(links = list(ndt = "log")), b),
               "takes the family object")
  expect_error(ndt_bound_attach(seam_family(), list(ub = 0.3)),
               "takes what ndt_bound\\(\\) returns")
  no_ndt <- frmtmb::frmtmb_family(
    "nondt", dpars = "mu", links = list(mu = "identity"),
    primary_dpars = "mu", type = "continuous",
    lpdf = function(y, dpars, aterms) 0)
  expect_error(ndt_bound_attach(no_ndt, b),
               "no distributional parameter named")
  # THIS PACKAGE'S OWN FAMILIES, ALL FIVE. The first version of this
  # guard tested `ndt_raw`, the marker ddm_ndt_install() leaves, and
  # only FOUR of the five carry it: gddm() passed, and the accepted call
  # replaced its bounded scaled logit with a plain logit, planted an
  # `ndt_floor` its density never reads, and moved its `ndt` starting
  # value from 0.155 s to 0.5. The guard now tests the property.
  for (fm in list(wiener(), lba(2), rdm(2), wiener_gng(), gddm())) {
    expect_error(ndt_bound_attach(fm, b), "own families",
                 info = fm[["family"]])
  }
  # and the case where the guarded thing is ABSENT: gddm() is the one
  # with no `ndt_raw`, so it is the one the old guard could not see
  expect_null(gddm()[["ndt_raw"]])
  expect_true("ndt" %in% names(gddm()[["links"]]))
  # what the old guard would have let through, stated as the damage
  expect_match(gddm()[["links"]][["ndt"]][["name"]], "scaled_logit")
})

test_that("ndt_bound refuses a response no bound can come from", {
  # ?ndt_bound invites a caller to assemble aterms by hand and call this
  # outside frm(), where no valid_y() runs in front of it. Each of these
  # came back with a bound before the checks were added, and the NA one
  # came back with the PENDING sentinel from a bound that was set.
  rt <- c(0.31, 0.42, 0.55)
  expect_error(ndt_bound("a", what = "f"), "positive, finite times")
  expect_error(ndt_bound(c(0.3, NA), what = "f"), "missing value")
  expect_error(ndt_bound(c(-0.3, 0.4), what = "f"), "at or below zero")
  expect_error(ndt_bound(c(0, 0.4), what = "f"), "at or below zero")
  expect_error(ndt_bound(numeric(0), what = "f"), "positive, finite")
  expect_error(ndt_bound(c(0.3, Inf), what = "f"), "positive, finite")
  # and the two entry points now agree about max_ndt, where they did not
  for (bad in list(0, -1, c(1, 2), NA_real_, "x")) {
    expect_error(ndt_bound(rt, max_ndt = bad, what = "f"),
                 "one positive number")
    expect_error(ndt_bound_pending(bad, what = "f"),
                 "one positive number")
  }
  # what stays accepted
  expect_identical(ndt_bound(rt, max_ndt = 0.2, what = "f")[["ub"]], 0.2)
  expect_identical(ndt_bound(rt, what = "f")[["ub"]], 0.31)
})

test_that("ndt_apply multiplies the floor out, or refuses", {
  # The sixth export, and the reason it is one: the arithmetic is three
  # lines and the refusal is longer, so a consumer who copies the first
  # and not the second reads a fraction as seconds.
  expect_identical(ndt_apply(list(ndt = 0.2)), 0.2)
  expect_equal(ndt_apply(list(ndt = 0.5), list(ndt_floor = c(0.4, 0.6))),
               c(0.2, 0.3))
  # one list serving as both, which is what a merged engine hands it
  expect_equal(ndt_apply(list(ndt = 0.5, ndt_floor = c(0.4, 0.6))),
               c(0.2, 0.3))
  # THE CASE WHERE THE GUARDED THING IS ABSENT: grouped, floor missing
  expect_error(ndt_apply(list(ndt = 0.5, ndt_group = c(1, 2)),
                         what = "myfam()"),
               "myfam[(][)]: this model bounds")
  # the grouping present WITH the floor is not the refusing case
  expect_equal(ndt_apply(list(ndt = 0.5, ndt_group = c(1, 2),
                              ndt_floor = c(0.4, 0.6))),
               c(0.2, 0.3))
  # and it tapes: no comparison and no branch on a parameter
  tp <- RTMB::MakeTape(function(p) {
    sum(ndt_apply(list(ndt = p), list(ndt_floor = c(0.4, 0.6))))
  }, 0.5)
  expect_equal(as.numeric(tp$jacobian(0.5)), 1.0)
})

test_that("a family that uses the seam end to end scores", {
  # the seam's own worked example, assembled here rather than described:
  # derive, attach, and let the density call ndt_apply().
  rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
  at <- list(ndt_group = ndt_bound_key(c("a", "a", "a", "b", "b")))
  fam <- ndt_bound_attach(seam_family(), ndt_bound(rt, at, what = "s"))
  av <- c(at, fam[["aterm_data"]](rt, at))
  # `ndt` is a fraction and the density gets a time
  got <- fam[["lpdf"]](rt, list(mu = 0, ndt = 0.5), av)
  want <- stats::dnorm(rt - c(0.155, 0.155, 0.155, 0.305, 0.305), 0, 1,
                       log = TRUE)
  expect_equal(as.numeric(got), want)
  # and with the floor gone it refuses instead of scoring the fraction
  expect_error(fam[["lpdf"]](rt, list(mu = 0, ndt = 0.5), at),
               "this model bounds")
})

test_that("the pending state refuses rather than reporting a wrong scale", {
  # bf(ndt = ) runs linkfun() at PARSE time and mixture() never
  # finalizes its components, so a family with an unsettled bound must
  # refuse rather than transform a number on a scale nothing has set.
  out <- ndt_bound_attach(seam_family(),
                          ndt_bound_pending(what = "seamtest"))
  expect_null(ndt_bound_of(out))
  expect_error(out$links$ndt$linkinv(0), "bound is not set yet")
  expect_error(out$links$ndt$linkfun(0.2), "bound is not set yet")
  settled <- ndt_bound_attach(seam_family(),
                              ndt_bound_pending(0.4, what = "seamtest"))
  expect_equal(settled$links$ndt$linkinv(0), 0.2)
  expect_identical(ndt_bound_of(settled)$ub, 0.4)
})

test_that("a family carrying the seam's bound is reported by ndt_time", {
  # the whole reason ndt_time() reads the RECORD rather than the
  # family's name: a fit from another package has to be reportable
  skip_if_not_installed("RWiener")
  set.seed(4242)
  d <- ddm_simulate(200, mu = 1.2, bs = 1.5, ndt = 0.25)
  d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))
  fit <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
                bias = 0.5),
             family = wiener(), data = d)
  bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  expect_s3_class(bd, "frmtmb_eam_ndt_bound")
  expect_identical(ndt_bound_of(frmtmb::single_response(fit)[["family"]]),
                   bd)
  expect_identical(bd$what, "wiener")
})
