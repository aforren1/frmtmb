## The family object, its components, and every refusal. A family is not
## finished when it fits; it is finished when the things a user gets wrong
## are refused in words that say what to do instead.

test_that("the component catalogue declares the parameters it says it does", {
  expect_setequal(gddm()[["dpars"]], c("mu", "bs", "bias", "ndt"))
  expect_identical(gddm()[["dpars"]][[1L]], "mu")

  f <- gddm(drift = list(gddm_drift_coherence(cmax = 0.512),
                         gddm_drift_leak()),
            bound = gddm_bound_exponential(),
            start = gddm_start_uniform(),
            lapse = "uniform")
  expect_setequal(f[["dpars"]],
                  c("mu", "alpha", "leak", "bs", "tau", "bias", "sz",
                    "ndt", "lapse"))
  expect_identical(f[["dpars"]][[1L]], "mu")
  expect_identical(f[["links"]][["alpha"]]$name, "log")
  expect_identical(f[["links"]][["leak"]]$name, "identity")
  expect_identical(f[["links"]][["bs"]]$name, "log")
  expect_identical(f[["links"]][["kappa"]], NULL)

  fl <- gddm(bound = gddm_bound_linear())
  expect_true("kappa" %in% fl[["dpars"]])
  expect_identical(fl[["links"]][["kappa"]]$name, "logit")
})

test_that("the family declares the addition terms it cannot do without", {
  ## Only the covariates can be declared. The boundary arrives as dec()
  ## OR vint1, and the condition's slot moves with that choice, so
  ## neither can be named in a declaration that means "all of these".
  expect_identical(gddm()[["required_aterms"]], character(0))
  expect_setequal(gddm(drift = gddm_drift_coherence())[["required_aterms"]],
                  "vreal1")
})

test_that("the boundary and the condition are read under either spelling", {
  ## dec() is the spelling brms uses and the one this package registers;
  ## vint() is the general-purpose route that predates the registry.
  ## vint() numbers positionally, so the condition moves when dec()
  ## takes over the boundary, and the reader has to move with it.
  viaint <- gd_indicator(list(vint1 = c(0, 1), vint2 = c(1L, 2L)))
  expect_identical(viaint[["up"]], c(0, 1))
  expect_identical(viaint[["cond"]], c(1L, 2L))

  viadec <- gd_indicator(list(dec = c(0, 1), vint1 = c(1L, 2L)))
  expect_identical(viadec[["up"]], c(0, 1))
  expect_identical(viadec[["cond"]], c(1L, 2L))

  ## dec() wins when both are present, and the message names where the
  ## condition then lives
  both <- gd_indicator(list(dec = c(1, 0), vint1 = c(3L, 4L),
                            vint2 = c(9L, 9L)))
  expect_identical(both[["up"]], c(1, 0))
  expect_identical(both[["cond"]], c(3L, 4L))
  expect_match(both[["cond_is"]], "first value of vint")
  expect_match(viaint[["cond_is"]], "second value of vint")

  ## and nothing at all reads as absent rather than as an error
  none <- gd_indicator(list())
  expect_null(none[["up"]])
  expect_null(none[["cond"]])
})

test_that("gddm() refuses a component in the wrong slot", {
  expect_error(gddm(drift = gddm_bound_constant()),
               "every element of `drift` must be a drift component")
  expect_error(gddm(bound = gddm_drift_constant()),
               "`bound` must be one boundary component")
  expect_error(gddm(start = gddm_bound_constant()),
               "`start` must be one starting-point component")
  expect_error(gddm(control = list(dt = 1)), "must come from gddm_control")
})

test_that("a drift needs exactly one base term, first", {
  expect_error(gddm(drift = gddm_drift_leak()), "needs exactly one base term")
  expect_error(gddm(drift = list(gddm_drift_leak(), gddm_drift_constant())),
               "needs exactly one base term")
  expect_error(gddm(drift = list(gddm_drift_constant(),
                                 gddm_drift_coherence())),
               "needs exactly one base term")
  expect_s3_class(gddm(drift = list(gddm_drift_constant(),
                                    gddm_drift_leak())),
                  "frmtmb_family")
})

test_that("gddm_control validates its own arguments", {
  expect_error(gddm_control(dt = 0), "`dt` must be one positive")
  expect_error(gddm_control(dt = c(0.1, 0.2)), "`dt` must be one positive")
  expect_error(gddm_control(ny = 3), "`ny` must be one integer of at least")
  expect_error(gddm_control(t_max = -1), "`t_max` must be one positive")
  expect_error(gddm_control(renormalize = NA), "must be TRUE or FALSE")
  expect_error(gddm_control(max_ndt = 0), "`max_ndt` must be one positive")
  expect_error(gddm_control(tridiagonal = "magic"), "should be one of")
  expect_s3_class(gddm_control(), "gddm_control")
  expect_identical(gddm_control()$tridiagonal, "recorded")
})

test_that("gddm_drift_coherence validates its scale", {
  expect_error(gddm_drift_coherence(cmax = 0), "`cmax` must be one positive")
  expect_error(gddm_drift_coherence(cov = 0), "`cov` must be one positive")
})

gd_toy <- function(n = 120, seed = 5) {
  set.seed(seed)
  d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                     control = gddm_control(t_max = 2))
  d$cond <- 1L
  d
}
gd_small <- function() gddm_control(t_max = 2, dt = 0.05, ny = 51L)

test_that("a decision indicator with more than two values is refused by name", {
  ## One accumulator between two absorbing boundaries can end a trial at
  ## the upper wall or the lower wall. A third response is a different
  ## model, not another parameter, and folding it into one of the two
  ## would silently fit something nobody asked for.
  d <- gd_toy()
  ## set outright rather than perturbing a draw, so the test pins three
  ## levels however the simulation happens to fall
  d$upper <- rep(c(0L, 1L, 2L), length.out = nrow(d))
  expect_error(
    frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
        family = gddm(control = gd_small()), data = d),
    "3 distinct values, and this family admits exactly two")
  expect_error(
    frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
        family = gddm(control = gd_small()), data = d),
    "racing accumulators")
  ## and it names the family that fits them, which lands in this same
  ## package: a refusal that says "not here" would be false after merge
  expect_error(
    frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
        family = gddm(control = gd_small()), data = d),
    "lba(n)", fixed = TRUE)

  ## four levels names four
  d$upper <- rep(c(0L, 1L, 2L, 3L), length.out = nrow(d))
  expect_error(
    frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
        family = gddm(control = gd_small()), data = d),
    "4 distinct values")
})

test_that("a two-valued indicator that is not 0/1 is refused with the recoding", {
  d <- gd_toy()
  d$upper <- rep(c(1L, 2L), length.out = nrow(d))
  expect_error(
    frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
        family = gddm(control = gd_small()), data = d),
    "must be 0 at the lower boundary and 1 at the upper one")
})

test_that("a missing addition term is refused, not silently defaulted", {
  d <- gd_toy()
  expect_error(frm(bf(rt ~ 1, bias = 0.5),
                   family = gddm(control = gd_small()), data = d),
               "decision indicator is missing")
  ## and it names both spellings, because either one supplies it
  expect_error(frm(bf(rt ~ 1, bias = 0.5),
                   family = gddm(control = gd_small()), data = d),
               "dec(", fixed = TRUE)
  expect_error(frm(bf(rt ~ 1, bias = 0.5),
                   family = gddm(control = gd_small()), data = d),
               "vint(upper, cond)", fixed = TRUE)
  ## a boundary with no condition beside it is refused on its own terms
  expect_error(frm(bf(rt | vint(upper) ~ 1, bias = 0.5),
                   family = gddm(control = gd_small()), data = d),
               "condition index is missing")
  ## the coherence drift needs its covariate too
  expect_error(
    frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
        family = gddm(drift = gddm_drift_coherence(cmax = 0.5),
                      control = gd_small()), data = d),
    "vreal")
})

test_that("a response time that is not a response time is refused", {
  d <- gd_toy()
  d$rt[3] <- -1
  expect_error(frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
                   family = gddm(control = gd_small()), data = d),
               "strictly positive, finite response time")
})

test_that("a covariate that is not constant within a condition is refused", {
  ## One solve serves a whole condition, so a covariate the drift reads
  ## has to take one value there. This is the part of the
  ## constant-within-condition contract the family can check.
  set.seed(6)
  d <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25, coh = c(0.1, 0.5),
                     drift = gddm_drift_coherence(cmax = 0.5),
                     control = gddm_control(t_max = 2))
  d$cond <- 1L
  expect_error(
    frm(bf(rt | vint(upper, cond) + vreal(coh) ~ 1, bias = 0.5),
        family = gddm(drift = gddm_drift_coherence(cmax = 0.5),
                      control = gd_small()), data = d),
    "not constant within every condition")
})

test_that("a window that does not contain the data is refused", {
  d <- gd_toy()
  expect_error(
    frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
        family = gddm(control = gddm_control(t_max = 0.2, dt = 0.05,
                                             ny = 51L)), data = d),
    "at or below the largest response time")
})

test_that("a non-decision-time bound above the fastest response is refused", {
  d <- gd_toy()
  expect_error(
    frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
        family = gddm(control = gddm_control(t_max = 2, dt = 0.05, ny = 51L,
                                             max_ndt = 10)), data = d),
    "above the smallest response time")
})

test_that("the links refuse to be read before the data resolves them", {
  ## The non-decision-time bound comes from the response, so the family
  ## says so rather than naming a link it will not use.
  expect_error(gddm()[["links"]][["ndt"]]$linkinv(0),
               "not resolved yet")
  ## unless the bound was pinned up front
  f <- gddm(control = gddm_control(max_ndt = 0.4))
  expect_error(f[["links"]][["ndt"]]$linkinv(0), "not resolved yet")
})

test_that("gddm_conditions numbers the distinct combinations", {
  d <- data.frame(coh = c(0, 0, 0.5, 0.5), block = c(1, 2, 1, 2))
  expect_identical(gddm_conditions(d, coh), c(1L, 1L, 2L, 2L))
  expect_identical(gddm_conditions(d, ~ coh), c(1L, 1L, 2L, 2L))
  expect_identical(gddm_conditions(d, coh, block), c(1L, 2L, 3L, 4L))
  expect_error(gddm_conditions(d), "name at least one variable")
  expect_error(gddm_conditions(d, ~ nosuch), "which `data` does not have")
  expect_error(gddm_conditions(1:3, coh), "must be a data frame")
})

test_that("gddm_simulate refuses parameters the components do not have", {
  expect_error(gddm_simulate(10, mu = 1, tau = 2),
               "have no parameter tau")
  expect_error(gddm_simulate(10, mu = 1, ndt = 5,
                             control = gddm_control(t_max = 2)),
               "at or past the end of the simulated window")
})

test_that("gddm_simulate draws from the model it is given", {
  set.seed(7)
  ## a positive drift favors the upper boundary, a negative one the lower
  up <- gddm_simulate(800, mu = 2.5, bs = 2, ndt = 0.2,
                      control = gddm_control(t_max = 2))
  dn <- gddm_simulate(800, mu = -2.5, bs = 2, ndt = 0.2,
                      control = gddm_control(t_max = 2))
  expect_gt(mean(up$upper), 0.8)
  expect_lt(mean(dn$upper), 0.2)
  ## every response time is past the non-decision time and inside the window
  expect_true(all(up$rt > 0.2))
  expect_true(all(up$rt <= 2))
  ## a wider boundary separation is slower
  wide <- gddm_simulate(800, mu = 2.5, bs = 3, ndt = 0.2,
                        control = gddm_control(t_max = 3))
  expect_gt(mean(wide$rt), mean(up$rt))
  expect_named(up, c("rt", "upper", "cond"))
})

test_that("the family prints what it is", {
  expect_output(print(gddm_drift_coherence()), "coherence")
  expect_output(print(gddm_drift_coherence()), "mu, alpha")
  expect_output(print(gddm_bound_exponential()), "bound component")
  expect_output(print(gddm_control()), "renormalize = TRUE")
})

test_that("every message this family raises is its own", {
  ## A shared message template means one of the two callers is telling the
  ## user about the wrong thing. This is the ddm clone of frmtmb's own
  ## message-uniqueness check, restricted to the generalized family.
  ##
  ## Read from the namespace, not from R/gddm.R: an installed package has
  ## no source tree, and this has to run under R CMD check as well as from
  ## the working copy.
  ns <- asNamespace("frmtmb.ddm")
  txt <- paste(vapply(ls(ns, all.names = TRUE), function(nm) {
    o <- get(nm, envir = ns)
    if (is.function(o)) paste(deparse(o), collapse = "\n") else ""
  }, character(1L)), collapse = "\n")
  ## the first literal of a call is its template; the rest of the call is
  ## the part that varies with the data
  starts <- regmatches(
    txt,
    gregexpr('(stop|warning|message)\\(\\s*"gddm[^"]{15,}"', txt))[[1L]]
  starts <- sub('^(stop|warning|message)\\(\\s*"', "", starts)
  expect_gt(length(starts), 8L)
  expect_identical(anyDuplicated(starts), 0L)
})
