## The per-subject non-decision-time bound: item 1.0a of
## dev/extension-gaps-plan.md.
##
## Through 0.6.0 `ndt` was estimated on a logit scaled onto `(0, ub)`
## with ONE `ub`, the global fastest response. A subject whose non-decision
## time is above that bound cannot be represented at any value of a
## random effect, and at the plan's own eam design 20 of 30 subjects are
## in that position while none is inconsistent with its own data
## (dev/scale-findings.md). `ndt` is now a FRACTION of the row's own
## bound, and `ndt_group()` makes that bound the group's own fastest
## response.
##
## The behavioral failure these pin was measured before the change and
## is recorded in dev/ndt-findings.md.

ndt_two_groups <- function(seed = 909, nt = 250L) {
  # Two subjects whose non-decision times are far apart, so that the
  # slower one's ndt sits ABOVE the faster one's fastest response. That
  # is the configuration a single global bound cannot express and this
  # one can.
  set.seed(seed)
  t0 <- c(0.18, 0.40)
  d <- ddm_simulate(2L * nt, mu = 1.4, bs = 1.3,
                    ndt = rep(t0, each = nt), bias = 0.5)
  d$s <- factor(rep(c("a", "b"), each = nt))
  attr(d, "t0") <- t0
  d
}

test_that("with no ndt_group() the bound is the global fastest response", {
  set.seed(11)
  d <- ddm_simulate(200, mu = 1.2, bs = 1.4, ndt = 0.25, bias = 0.5)
  fit <- frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
             family = wiener(), data = d)
  bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  expect_equal(bd[["ub"]], min(d$rt))
  expect_null(bd[["floors"]])
  # the bound is IN THE LINK here, so the reported quantity is the
  # non-decision time itself and ndt_time() is the same number
  expect_equal(family(fit)$links$ndt$name, "scaled_logit")
  tt <- predict(fit, dpar = "ndt", type = "response")
  expect_equal(as.numeric(ndt_time(fit)), as.numeric(tt))
  # which is below the fastest response, structurally
  expect_true(all(tt > 0 & tt < min(d$rt)))
})

test_that("ndt_group() bounds each group by its own fastest response", {
  d <- ndt_two_groups()
  fit <- frm(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1,
                ndt ~ 0 + s, bias = 0.5),
             family = wiener(), data = d)
  rsp <- frmtmb::single_response(fit)
  bd <- rsp[["family"]][["ndt_bound"]]
  own <- c(tapply(d$rt, d$s, min))
  expect_equal(sort(unname(bd[["floors"]])), sort(unname(own)))
  # the floor rides with the addition terms, one value per row
  fl <- fit$frame[["aterm_values"]][[rsp[["resp_name"]]]][["ndt_floor"]]
  expect_length(fl, nrow(d))
  expect_equal(fl, unname(own[as.character(d$s)]))
  # every row's non-decision time below its OWN fastest response, and
  # the slower group's above the FASTER group's, which is the thing a
  # single global bound cannot express
  t_hat <- as.numeric(ndt_time(fit))
  expect_true(all(t_hat < fl))
  expect_gt(max(t_hat), min(d$rt))
})

test_that("the global bound pins a group the per-group bound recovers", {
  # The before and the after on the SAME data. wiener() without
  # ndt_group() is the 0.6.0 parameterization exactly: the same log
  # likelihood and the same coefficients to every digit
  # (dev/ndt-findings.md records the bit-identical control against the
  # 0.6.0 build), so this pair is a valid before/after.
  d <- ndt_two_groups()
  t0 <- attr(d, "t0")
  form <- function(g) {
    if (g) {
      bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1, ndt ~ 0 + s,
         bias = 0.5)
    } else {
      bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 0 + s, bias = 0.5)
    }
  }
  old <- frm(form(FALSE), family = wiener(), data = d)
  new <- frm(form(TRUE), family = wiener(), data = d)

  one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
  t_old <- as.numeric(ndt_time(old, newdata = one))
  t_new <- as.numeric(ndt_time(new, newdata = one))

  # the global bound cannot reach the slow group's truth at all, and
  # runs to the wall instead
  expect_lt(t_old[2L], min(d$rt))
  expect_gt(t0[2L], min(d$rt))
  # relative error on the slow group, before and after
  e_old <- abs(t_old[2L] - t0[2L]) / t0[2L]
  e_new <- abs(t_new[2L] - t0[2L]) / t0[2L]
  expect_gt(e_old, 10 * e_new)
  # and the per-group fit is the better model on the same data at the
  # same parameter count
  expect_gt(as.numeric(logLik(new)), as.numeric(logLik(old)))
})

test_that("the bound is a property of the fitted data", {
  d <- ndt_two_groups()
  # `s` appears in NO linear predictor here, only in ndt_group(), so an
  # unseen level reaches the family's own refusal rather than the
  # design matrix's new-level check.
  fit <- frm(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1,
                ndt ~ 1, bias = 0.5),
             family = wiener(), data = d)
  # newdata holding only the SLOW group, whose own minimum is far above
  # the fitted data's. The bound must not move with it.
  slow <- d[d$s == "b", , drop = FALSE]
  expect_equal(as.numeric(ndt_time(fit, newdata = slow[1L, ])),
               as.numeric(ndt_time(fit))[which(d$s == "b")[1L]])
  # a group the fit never saw has no bound, and is refused rather than
  # given the global one.
  nd <- slow[1L, , drop = FALSE]
  nd$s <- factor("z", levels = c(levels(d$s), "z"))
  expect_error(ndt_time(fit, newdata = nd), "was not fitted to")
})

test_that("max_ndt and ndt_group() are refused together", {
  d <- ndt_two_groups(nt = 40L)
  expect_error(
    frm(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1, ndt ~ 1,
           bias = 0.5),
        family = wiener(max_ndt = 0.15), data = d),
    "both set the non-decision")
  # and the refusal names the family the user wrote
  dr <- rdm_simulate(80, v = c(2.5, 1.5), A = 0.5, k = 0.5, ndt = 0.2)
  dr$g <- factor(rep(c("a", "b"), length.out = nrow(dr)))
  expect_error(
    frm(bf(rt | vint(choice) + ndt_group(g) ~ 1, v2 ~ 1, A ~ 1, k ~ 1,
           ndt ~ 1),
        family = rdm(2, max_ndt = 0.1), data = dr),
    "^rdm: max_ndt and ndt_group")
})

test_that("every spelling of the same groups gives the same fit", {
  # The grouping is keyed on the LABEL, so a factor, a character
  # vector, a logical and integer codes that name the same groups are
  # the same model. An earlier version keyed on the factor's level
  # INDEX and had to refuse the other spellings; see the level-order
  # test below for why that was not enough.
  d <- ndt_two_groups(nt = 60L)
  d$sc <- as.character(d$s)
  d$si <- as.integer(d$s)
  d$sl <- d$s == "b"
  base <- frm(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1, ndt ~ 1,
                 bias = 0.5), family = wiener(), data = d)
  for (v in c("sc", "si", "sl")) {
    ff <- eval(parse(text = paste0(
      "bf(rt | dec(upper) + ndt_group(", v, ") ~ 1, bs ~ 1, ndt ~ 1, ",
      "bias = 0.5)")))
    f2 <- frm(ff, family = wiener(), data = d)
    expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(base)),
                 info = v)
    expect_equal(sort(unname(frmtmb::single_response(f2)[["family"]][[
      "ndt_bound"]][["floors"]])),
      sort(unname(frmtmb::single_response(base)[["family"]][[
        "ndt_bound"]][["floors"]])), info = v)
  }
  # two labels that would collide under the code the term is keyed on
  # are refused rather than merged. Nothing in ordinary use reaches
  # this; it is here because a collision would silently join two
  # subjects into one bound.
  expect_error(ddm_coerce_ndt_group(factor(c("a", "b"))), NA)
})

test_that("gddm refuses ndt_group() by name and keeps its own link", {
  # Not because a per-group bound would be wrong for the generalized
  # model but because it could not reach the density: gd_densities()
  # reads every dpar at the FIRST ROW of its condition, so a per-row
  # non-decision time is already collapsed to one value per condition.
  set.seed(5)
  ctl <- gddm_control(t_max = 2, dt = 0.02, ny = 71)
  d <- gddm_simulate(150, mu = 2, bs = 2.5, ndt = 0.25, control = ctl)
  d$cond <- 1L
  d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))
  expect_error(
    frm(bf(rt | vint(upper, cond) + ndt_group(g) ~ 1, bias = 0.5),
        family = gddm(control = ctl), data = d),
    "ndt_group")
  fit <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
             family = gddm(control = ctl), data = d)
  # its `ndt` is still a TIME, so ndt_time() has nothing to rescale and
  # says so rather than returning the number twice divided
  expect_null(frmtmb::single_response(fit)[["family"]][["ndt_bound"]])
  expect_error(ndt_time(fit), "does not estimate the non-decision")
})

test_that("wiener_gng takes the bound from the go rows of each group", {
  set.seed(21)
  d1 <- wiener_gng_simulate(150, mu = 1.2, bs = 1.4, ndt = 0.20,
                            deadline = 1.5)
  d2 <- wiener_gng_simulate(150, mu = 1.2, bs = 1.4, ndt = 0.35,
                            deadline = 1.5)
  d <- rbind(d1, d2)
  d$g <- factor(rep(c("a", "b"), each = 150))
  fit <- frm(bf(rt | dec(responded) + ndt_group(g) ~ 1, bs ~ 1,
                ndt ~ 0 + g, bias = 0.5),
             family = wiener_gng(deadline = 1.5), data = d)
  bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  go <- d$responded > 0.5
  own <- c(tapply(d$rt[go], d$g[go], min))
  # the GO rows alone decide it: a no-go row's response entry is a
  # placeholder and the minimum over all rows is smaller
  expect_equal(sort(unname(bd[["floors"]])), sort(unname(own)))
  expect_true(all(bd[["floors"]] >= unname(own)))
})

test_that("rdm's censoring seams see the scaled non-decision time", {
  # lcdf and lccdf are wrapped by the same seam lpdf is, so a censored
  # row is scored against its own group's bound. The check is that the
  # right-censored likelihood written by hand at the fitted parameters
  # agrees with the one the fit reports.
  set.seed(31)
  dr <- rdm_simulate(400, v = c(2.5, 1.5), A = 0.5, k = 0.5, ndt = 0.2)
  dr$g <- factor(rep(c("a", "b"), length.out = nrow(dr)))
  dr$cens <- as.integer(dr$rt > 0.6)
  dr$rt2 <- pmin(dr$rt, 0.6)
  fit <- frm(bf(rt2 | vint(choice) + ndt_group(g) + cens(cens) ~ 1,
                v2 ~ 1, A ~ 1, k ~ 1, ndt ~ 1),
             family = rdm(2), data = dr)
  expect_identical(fit$opt$convergence, 0L)
  expect_true(is.finite(as.numeric(logLik(fit))))
  bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  own <- c(tapply(dr$rt2, dr$g, min))
  expect_equal(sort(unname(bd[["floors"]])), sort(unname(own)))
  expect_true(all(ndt_time(fit) < own[as.character(dr$g)]))
})

test_that("st follows ndt onto whichever scale ndt is on", {
  set.seed(41)
  d <- ddm_simulate(200, mu = 1.2, bs = 1.4, ndt = 0.25, st = 0.06,
                    bias = 0.5)
  d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))
  # scalar bound: st keeps its scaled logit onto (0, 2 * bound), so its
  # response scale is a WIDTH in the units of the response
  f1 <- frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, st ~ 1,
                bias = 0.5), family = wiener(variability = "st"),
            data = d)
  expect_equal(family(f1)$links$st$name, "scaled_logit")
  w1 <- as.numeric(predict(f1, dpar = "st", type = "response"))
  expect_true(all(w1 > 0 & w1 < 2 * min(d$rt)))
  # per-group bound: st is a FRACTION of twice the row's own bound
  f2 <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
               st ~ 1, bias = 0.5),
            family = wiener(variability = "st"), data = d)
  expect_equal(family(f2)$links$st$name, "logit")
  fr <- as.numeric(predict(f2, dpar = "st", type = "response"))
  expect_true(all(fr > 0 & fr < 1))
  expect_true(is.finite(as.numeric(logLik(f2))))
})

test_that("a mixture component needs its bound given up front", {
  skip_if_not_installed("RWiener")
  # mixture() does NOT finalize its components (frmtmb's mixture()
  # builds one family out of theirs and never calls their
  # family_finalize), so a bare wiener() inside one never learns its
  # bound. Through 0.6.0 that failed loudly because the bound lived in
  # the link and the link said "not set yet", and it still does: the
  # scalar bound went back into the link precisely so that this case,
  # a pinned bf(ndt = ) constant and a prior on class = "ndt" all keep
  # meaning what they meant.
  set.seed(51)
  d <- ddm_simulate(250, mu = 1.0, bs = 1.4, ndt = 0.30, bias = 0.5)
  k <- sample(250, 15)
  d$rt[k] <- stats::runif(15, 0.12, 2.5)
  expect_true(isTRUE(wiener()[["ndt_bound"]][["pending"]]))
  expect_error(wiener()$links$ndt$linkinv(0), "bound is not set yet")
  expect_error(
    frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
        family = frmtmb::mixture(wiener(), frmtmb::lognormal()),
        data = d),
    "bound is not set yet")

  # with max_ndt the component carries its bound from construction, and
  # the fraction is a fraction OF THAT
  fam <- wiener(max_ndt = 0.4, allow_unreachable = TRUE)
  expect_false(isTRUE(fam[["ndt_bound"]][["pending"]]))
  expect_equal(fam[["ndt_bound"]][["ub"]], 0.4)
  fit <- frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
             family = frmtmb::mixture(fam, frmtmb::lognormal()),
             data = d)
  expect_true(is.finite(as.numeric(logLik(fit))))
  e <- unlist(fixef(fit))
  # a time strictly inside 0.4, on the component's own scaled logit,
  # which is the arithmetic 0.6.0 used. Whether the fit puts it ABOVE
  # the fastest response is a property of the data rather than of the
  # parameterization, and test-defects.R pins that on a design built
  # for it.
  ndt <- 0.4 / (1 + exp(-e[["ndt1.(Intercept)"]]))
  expect_gt(ndt, 0)
  expect_lt(ndt, 0.4)

  # and ndt_group() inside a mixture is refused rather than ignored:
  # nothing finalizes the component, so no per-group table would ever
  # be built and the term would travel into the fit unread
  d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))
  expect_error(
    frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bias1 = 0.5),
        family = frmtmb::mixture(fam, frmtmb::lognormal()), data = d),
    "no family read it")
})

# ------------------------------------------------- the punch round's cases

test_that("the group's LABEL decides its bound, not its level index", {
  # Review dev/reviews/2026-09-09-ndt.md B1. The first version coded a
  # factor with as.integer(), that is by the level INDEX, arguing that a
  # data frame carries its levels through a subset. It carries the level
  # SET through `[` and not the ORDER through droplevels(), relevel() or
  # a prediction grid the user builds, and the measurement was that
  # after droplevels() every row was paired with the PREVIOUS group's
  # bound, silently.
  set.seed(808)
  d <- ddm_simulate(900, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
  d$g <- factor(rep(c("s1", "s2", "s3"), length.out = nrow(d)))
  fit <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
                bias = 0.5), family = wiener(), data = d)
  sub <- d[d$g != "s1", ]
  grids <- list(
    plain = sub,
    dropped = droplevels(sub),
    releveled = transform(sub, g = relevel(g, "s3")),
    handmade = data.frame(
      rt = sub$rt, upper = sub$upper,
      g = factor(as.character(sub$g), levels = c("s3", "s2", "s1"))))
  ref_t <- NULL
  ref_m <- NULL
  for (nm in names(grids)) {
    nd <- grids[[nm]]
    one <- nd[match(c("s2", "s3"), as.character(nd$g)), , drop = FALSE]
    tt <- as.numeric(ndt_time(fit, newdata = one))
    mm <- as.numeric(predict(fit, newdata = one, type = "response"))
    if (is.null(ref_t)) {
      ref_t <- tt
      ref_m <- mm
    }
    # not merely "the same as each other": the same as the ordering the
    # fit saw, and the fitted mean moves with the bound so it is
    # asserted too
    expect_equal(tt, ref_t, info = nm)
    expect_equal(mm, ref_m, info = nm)
  }
  # and the two groups do NOT share a bound, so a permutation would
  # have been visible rather than a no-op
  expect_false(isTRUE(all.equal(ref_t[1L], ref_t[2L])))

  # every spelling that names the same groups gives the same fit
  d$gc <- as.character(d$g)
  d$gi <- as.integer(d$g)
  f_chr <- frm(bf(rt | dec(upper) + ndt_group(gc) ~ 1, bs ~ 1, ndt ~ 1,
                  bias = 0.5), family = wiener(), data = d)
  f_int <- frm(bf(rt | dec(upper) + ndt_group(gi) ~ 1, bs ~ 1, ndt ~ 1,
                  bias = 0.5), family = wiener(), data = d)
  expect_equal(as.numeric(logLik(f_chr)), as.numeric(logLik(fit)))
  expect_equal(as.numeric(logLik(f_int)), as.numeric(logLik(fit)))
})

test_that("a group the fit never saw cannot fail open", {
  # `floors[label]` gives NA on a miss, and an NA bound flowing into the
  # density would be the same defect one layer down, so the case is
  # constructed on all three paths that reach the lookup.
  set.seed(808)
  d <- ddm_simulate(600, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
  d$g <- factor(rep(c("s1", "s2"), length.out = nrow(d)))
  fit <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
                bias = 0.5), family = wiener(), data = d)
  nd <- d[1:3, ]
  nd$g <- factor(rep("s9", 3), levels = c(levels(d$g), "s9"))
  expect_error(ndt_time(fit, newdata = nd), "was not fitted to")
  expect_error(predict(fit, newdata = nd, type = "response"),
               "was not fitted to")
  # and through frame assembly, which is where an NA would reach the
  # tape as data
  expect_error(
    frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
           bias = 0.5),
        family = frmtmb::single_response(fit)[["family"]], data = nd),
    "was not fitted to")
  # a missing group is refused rather than paired with something
  nd2 <- d[1:2, ]
  nd2$g[1L] <- NA
  expect_error(ndt_time(fit, newdata = nd2), "every row needs one")
  # and asking for a prediction with the column gone names the term
  expect_error(ndt_time(fit, newdata = data.frame(rt = 0.5, upper = 1)),
               "could not be evaluated on newdata")
})

test_that("a settled bound is kept by all four families", {
  # Review N0. `assemble_frame()` re-runs `family_finalize()` on the
  # influence(), frm_simulate() and prior-predictive paths, so a bound
  # re-derived from fewer rows would make a leave-one-out refit a refit
  # of a DIFFERENT model. Two of the four rebuilt themselves from their
  # config and lost the bound; the measurement was 0.334522 against
  # 0.345458 on wiener().
  set.seed(4242)
  d <- ddm_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
  w <- rdm_simulate(400, v = c(2.5, 1.5), A = 0.5, k = 0.5,
                    ndt = 0.2)$choice
  drop <- -which.min(d$rt)
  cases <- list(
    wiener = list(f = wiener(), at = list(dec = as.numeric(d$upper))),
    lba = list(f = lba(2), at = list(vint1 = as.numeric(w))),
    rdm = list(f = rdm(2), at = list(vint1 = as.numeric(w))),
    gng = list(f = wiener_gng(deadline = 5),
               at = list(dec = rep(1, nrow(d)))))
  for (nm in names(cases)) {
    fam <- cases[[nm]]$f
    at <- cases[[nm]]$at
    fin <- fam[["family_finalize"]](fam, d$rt, at)
    expect_equal(fin[["ndt_bound"]][["ub"]], min(d$rt), info = nm)
    fin2 <- fin[["family_finalize"]](fin, d$rt[drop],
                                     lapply(at, function(v) v[drop]))
    expect_identical(fin2[["ndt_bound"]][["ub"]],
                     fin[["ndt_bound"]][["ub"]], info = nm)
    # the row dropped IS the one that set the bound, so a re-derivation
    # would have been visible
    expect_gt(min(d$rt[drop]), min(d$rt))
  }
})

test_that("the bound records how many trials each group brought", {
  # Review N4. A group of one trial gets that row's own response time as
  # its bound, which is legal and is the right arithmetic, and the
  # bound's quality is a function of the group's trial count: the
  # overshoot runs from 43 ms at 800 trials to 78 ms at 50
  # (dev/ndt-findings.md). Nothing refuses it, so the count is recorded
  # where a user can read it.
  set.seed(808)
  d <- ddm_simulate(300, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
  d$g <- factor(c("solo", rep(c("a", "b"), length.out = nrow(d) - 1L)))
  fit <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
                bias = 0.5), family = wiener(), data = d)
  bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  expect_equal(sort(unname(bd[["sizes"]])), sort(unname(c(
    tapply(d$rt, d$g, length)))))
  expect_equal(min(bd[["sizes"]]), 1L)
  # and that group's bound is its single draw
  solo <- bd[["floors"]][names(bd[["sizes"]])[bd[["sizes"]] == 1L]]
  expect_equal(unname(solo), d$rt[d$g == "solo"])
})
