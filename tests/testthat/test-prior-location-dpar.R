# Where a prior names no dpar, no resp or a slot the model lacks, against
# brms 2.23.0's stancode() on the same models
# (dev/mvprior-log/brms-probe.txt, script dev/mvprior-brms-probe.R, data
# dev/mvprior-cases.R, seed 2309):
#
#   bf(cat ~ x) + categorical(): b, b coef x, Intercept: refused;
#     b dpar = "mub": lprior += normal_lpdf(b_mub | 0, 5)
#   bf(ym ~ x) + mixture(gaussian, gaussian): the same with mu1, mu2
#   bf(cat ~ x + (1 | g)): sd and sd group g refused; sd dpar = "mub"
#     reaches sd_1; default_prior() lists b, Intercept and sd per dpar
#   bf(y1 ~ x) + bf(cat ~ x, categorical()): b resp = "cat" refused
#   y ~ 1: class b refused; bf(y1 ~ x) + bf(y2 ~ 1): b resp = "y2" refused
#   a univariate model: b, sd and ar with resp = "y" refused
#   set_rescor(TRUE): rescor with resp refused
#   a ~ 1 + z in a nonlinear model: class b with no nlpar refused
#
# Through 0.62.0 frmtmb broadcast the first four to every location
# parameter, applied class b on y ~ 1 to nothing, applied rescor with a
# resp to the whole matrix, and put class b with no nlpar on a_z.
# The last three decisions are the user's, 2026-09-24.

ld_data <- function() {
  set.seed(2309)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(1:8, length.out = n)),
                  t = rep(1:10, each = 8))
  d$y <- 1 + d$x + stats::rnorm(8)[d$g] + stats::rnorm(n)
  d$y1 <- d$y
  d$y2 <- -1 + 0.5 * d$x + stats::rnorm(n)
  d$cat <- factor(c("a", "b", "c")[1 + (d$x + stats::rlogis(n) > 0) +
                                      (d$z + stats::rlogis(n) > 0.5)])
  d$ym <- ifelse(stats::rbinom(n, 1, 0.5) == 1, 3 + d$x, -1 + d$x) +
    stats::rnorm(n, 0, 0.5)
  d$ax <- abs(d$x)
  d$ynl <- (2 + 0.3 * d$z) * exp(-0.5 * d$ax) + stats::rnorm(n, 0, 0.1)
  d
}

ld_reached <- function(formula, data, pl) {
  des <- frmtmb:::prior_design(formula, data, NULL, list())
  r <- frmtmb:::resolve_priorlist(des, pl)
  pt <- des$frame[["par_template"]]
  nm <- unlist(lapply(r$entries, function(e) {
    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
  }))
  sort(c(nm, names(r$lower), names(r$upper)))
}

# the message of a refusal, or "" when there was none, so a block can
# go on asserting after a call that 0.62.0 accepts
ld_msg <- function(expr) {
  e <- tryCatch({
    force(expr)
    NULL
  }, error = identity)
  if (inherits(e, "frmtmb_error")) conditionMessage(e) else ""
}

several <- "location is several distributional parameters"

test_that("categorical: b and Intercept with no dpar are refused", {
  d <- ld_data()
  f <- bf(cat ~ x) + categorical()
  for (p in list(set_prior("normal(0, 5)", class = "b"),
                 set_prior("normal(0, 5)", class = "b", coef = "x"),
                 set_prior("normal(0, 5)", class = "Intercept"))) {
    expect_match(ld_msg(validate_prior(p, f, data = d)), several)
  }
  expect_match(ld_msg(validate_prior(set_prior("normal(0, 5)",
                                               class = "b"), f,
                                     data = d)), "mub, muc", fixed = TRUE)
})

test_that("mixture: b and Intercept with no dpar are refused", {
  d <- ld_data()
  f <- bf(ym ~ x) + mixture(gaussian(), gaussian())
  for (p in list(set_prior("normal(0, 5)", class = "b"),
                 set_prior("normal(0, 5)", class = "Intercept"))) {
    expect_match(ld_msg(validate_prior(p, f, data = d)), several)
  }
})

test_that("several locations: sd with no dpar is refused", {
  d <- ld_data()
  for (f in list(bf(cat ~ x + (1 | g)) + categorical(),
                 bf(ym ~ x + (1 | g)) + mixture(gaussian(), gaussian()))) {
    expect_match(ld_msg(validate_prior(set_prior("normal(0, 5)",
                                                 class = "sd"), f,
                                       data = d)),
                 "No random-effect SDs")
  }
})

test_that("several locations: a prior with dpar reaches that one", {
  d <- ld_data()
  f <- bf(cat ~ x + (1 | g)) + categorical()
  expect_identical(ld_reached(f, d, set_prior("normal(0, 5)", class = "b",
                                              dpar = "mub")), "mub_x")
  fm <- bf(ym ~ x) + mixture(gaussian(), gaussian())
  expect_identical(ld_reached(fm, d, set_prior("normal(0, 5)",
                                               class = "Intercept",
                                               dpar = "mu2")),
                   "mu2_(Intercept)")
  # last: 0.62.0 refuses it, and an error ends the block
  expect_length(ld_reached(f, d, set_prior("normal(0, 5)", class = "sd",
                                           dpar = "muc")), 1L)
})

test_that("categorical ~ 1 takes Intercept with dpar", {
  d <- ld_data()
  expect_identical(ld_reached(bf(cat ~ 1) + categorical(), d,
                              set_prior("normal(0, 5)", class = "Intercept",
                                        dpar = "mub")),
                   "mub_(Intercept)")
})

test_that("default_prior() lists several locations by dpar, as brms does", {
  d <- ld_data()
  tab <- as.data.frame(default_prior(bf(cat ~ x + (1 | g)) + categorical(),
                                     data = d))
  lab <- paste(tab$class, tab$coef, tab$group, tab$dpar, sep = "|")
  expect_true(all(c("b|||mub", "b|x||mub", "Intercept|||muc", "sd|||mub",
                    "sd||g|muc") %in% lab))
  expect_false(any(tab$class %in% c("b", "Intercept", "sd") &
                     !nzchar(tab$dpar)))
  tm <- as.data.frame(default_prior(bf(y1 ~ x, family = gaussian()) +
                                      bf(cat ~ x, family = categorical()) +
                                      set_rescor(FALSE), data = d))
  expect_true(any(tm$class == "b" & tm$resp == "cat" & tm$dpar == "mub"))
  expect_false(any(tm$class == "b" & tm$resp == "cat" & !nzchar(tm$dpar)))
})

test_that("a multivariate categorical response needs dpar and resp", {
  d <- ld_data()
  f <- bf(y1 ~ x, family = gaussian()) +
    bf(cat ~ x, family = categorical()) + set_rescor(FALSE)
  expect_match(ld_msg(validate_prior(set_prior("normal(0, 5)", class = "b",
                                               resp = "cat"), f, data = d)),
               several)
  expect_identical(ld_reached(f, d, set_prior("normal(0, 5)", class = "b",
                                              dpar = "mub", resp = "cat")),
                   "cat_mub_x")
})

test_that("the resp advice names every call that resolves, and only those", {
  d <- ld_data()
  f <- bf(y1 ~ x, family = gaussian()) +
    bf(cat ~ x, family = categorical()) + set_rescor(FALSE)
  m <- ld_msg(validate_prior(set_prior("normal(0, 5)", class = "b",
                                       dpar = "mub"), f, data = d))
  expect_match(m, "dpar = \"mub\", resp = \"cat\")", fixed = TRUE)
  expect_no_match(m, "No response of this model has this slot",
                  fixed = TRUE)
  fnl <- bf(y1 ~ x) +
    bf(ynl ~ a * exp(-b * ax), a ~ 1 + z, b ~ 1, nl = TRUE) +
    set_rescor(FALSE) + gaussian()
  m <- ld_msg(validate_prior(set_prior("normal(0, 5)", class = "Intercept",
                                       nlpar = "a"), fnl, data = d))
  expect_match(m, "nlpar = \"a\", resp = \"ynl\")", fixed = TRUE)
  expect_identical(ld_reached(fnl, d, set_prior("normal(0, 5)",
                                                class = "Intercept",
                                                nlpar = "a", resp = "ynl")),
                   "ynl_a_(Intercept)")
})

test_that("class b with no slope to reach is refused", {
  d <- ld_data()
  for (p in list(set_prior("normal(0, 1)", class = "b"),
                 set_prior("", class = "b", lb = 0))) {
    expect_match(ld_msg(validate_prior(p, bf(y ~ 1 + (1 | g)), data = d)),
                 "no population-level slope")
  }
  f <- bf(y1 ~ x) + bf(y2 ~ 1) + set_rescor(FALSE)
  expect_match(ld_msg(validate_prior(set_prior("normal(0, 1)", class = "b",
                                               resp = "y2"), f, data = d)),
               "no population-level slope")
  # the fit refuses it before fitting, so prior_summary() never lists it
  expect_match(ld_msg(frm(bf(y ~ 1) + gaussian(), data = d,
                          prior = set_prior("normal(0, 1)", class = "b"))),
               "no population-level slope")
})

test_that("resp on a univariate model is refused", {
  d <- ld_data()
  for (p in list(set_prior("normal(0, 5)", class = "b", resp = "y"),
                 set_prior("normal(0, 5)", class = "sd", resp = "y"),
                 set_prior("student_t(3, 0, 2.5)", class = "sigma",
                           resp = "y"))) {
    expect_match(ld_msg(validate_prior(p, bf(y ~ x + (1 | g)), data = d)),
                 "applies only to a multivariate model")
  }
  far <- bf(y ~ x + ar(time = t, gr = g, cov = TRUE))
  expect_match(ld_msg(validate_prior(set_prior("normal(0, 0.5)",
                                               class = "ar", resp = "y"),
                                     far, data = d)),
               "applies only to a multivariate model")
})

test_that("rescor with resp is refused", {
  d <- ld_data()
  f <- bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)
  expect_match(ld_msg(validate_prior(set_prior("lkj(2)", class = "rescor",
                                               resp = "y1"), f, data = d)),
               "takes no resp")
})

test_that("nonlinear: class b with no nlpar is refused, slopes or not", {
  d <- ld_data()
  f <- bf(ynl ~ a * exp(-b * ax), a ~ 1 + z, b ~ 1, nl = TRUE)
  for (p in list(set_prior("normal(0, 5)", class = "b"),
                 set_prior("", class = "b", lb = 0))) {
    expect_match(ld_msg(validate_prior(p, f, data = d)),
                 "location is nonlinear")
  }
})

# frmtmb's own spellings, kept by the user's decision of 2026-09-24
# although brms refuses both. They pass on 0.62.0 too.
test_that("cor with resp and Intercept with nlpar stay accepted", {
  d <- ld_data()
  f <- bf(mvbind(y1, y2) ~ x + (1 + x | g)) + set_rescor(FALSE)
  expect_length(ld_reached(f, d, set_prior("lkj(2)", class = "cor",
                                           resp = "y1")), 1L)
  fnl <- bf(ynl ~ a * exp(-b * ax), a ~ 1 + z, b ~ 1, nl = TRUE)
  expect_identical(ld_reached(fnl, d, set_prior("normal(0, 5)",
                                                class = "Intercept",
                                                nlpar = "a")),
                   "a_(Intercept)")
})

test_that("an unknown resp is named as such, not as a missing slope", {
  d <- ld_data()
  m <- ld_msg(validate_prior(set_prior("normal(0, 1)", class = "b",
                                       resp = "zzz"),
                             bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(FALSE),
                             data = d))
  expect_match(m, "names no response of this model. It has y1, y2",
               fixed = TRUE)
})

# cs() on the sequential ordinal families: brms 2.23.0 puts class "b" on
# the category-specific coefficients too, and names each term's row by
# its variable (dev/mvprior-log/cs-brms.txt, script dev/mvprior-cs.R,
# seed 4417):
#   acat ord ~ cs(x), class b:      normal_lpdf(to_vector(bcs) | 0, 0.123)
#   acat ord ~ z + cs(x), class b:  normal_lpdf(b ...) and to_vector(bcs)
#   coef = "x":                     normal_lpdf(bcs[1] | 0, 0.123)
#   the same for sratio and cratio; default_prior() lists b coef x
# Through 0.62.0 class b missed them: nothing on ord ~ cs(x), z alone on
# ord ~ z + cs(x), and coef = "x" was refused.
cs_data <- function() {
  set.seed(4417)
  n <- 150
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  w = stats::rnorm(n))
  d$ord <- factor(cut(d$x + 0.5 * d$z + stats::rlogis(n),
                      c(-Inf, -1, 0, 1.2, Inf), labels = FALSE),
                  ordered = TRUE)
  d
}

test_that("class b reaches the cs() coefficients, as in brms", {
  d <- cs_data()
  for (fam in list(acat(), sratio(), cratio())) {
    expect_identical(ld_reached(bf(ord ~ cs(x)) + fam, d,
                                set_prior("normal(0, 1)", class = "b")),
                     paste0("bcs2_", 1:3))
    expect_identical(ld_reached(bf(ord ~ z + cs(x)) + fam, d,
                                set_prior("normal(0, 1)", class = "b")),
                     c(paste0("bcs2_", 1:3), "z"))
  }
})

# one assertion per block: 0.62.0 refuses both, and an error ends a block
cs_coef_reached <- function(coef) {
  dd <- frmtmb:::prior_design(ord ~ z + cs(x) + cs(w), cs_data(), acat(),
                              list())
  r <- frmtmb:::resolve_priorlist(dd, set_prior("normal(0, 1)", class = "b",
                                                coef = coef))
  unlist(lapply(r$entries, function(e) paste0(e$comp, "_", e$idx)))
}

test_that("coef names the second cs() term, as in brms", {
  expect_setequal(cs_coef_reached("w"), paste0("bcs3_", 1:3))
})

test_that("coef names the first cs() term, as in brms", {
  expect_setequal(cs_coef_reached("x"), paste0("bcs2_", 1:3))
})

test_that("default_prior() lists a cs() term under class b", {
  d <- cs_data()
  tab <- as.data.frame(default_prior(ord ~ cs(x), data = d,
                                     family = sratio()))
  expect_true(any(tab$class == "b" & tab$coef == "x"))
  expect_true(any(tab$class == "b" & !nzchar(tab$coef)))
})

test_that("one location dpar is named as one", {
  d <- ld_data()
  d$cat2 <- factor(ifelse(d$x > 0, "yes", "no"))
  m <- ld_msg(validate_prior(set_prior("normal(0, 1)", class = "b"),
                             bf(cat2 ~ x) + categorical(), data = d))
  expect_match(m, "the distributional parameter muyes", fixed = TRUE)
  expect_no_match(m, "several", fixed = TRUE)
})

test_that("the resp advice offers resp and dpar for a categorical response", {
  d <- ld_data()
  f <- bf(y1 ~ x, family = gaussian()) +
    bf(cat ~ x, family = categorical()) + set_rescor(FALSE)
  m <- ld_msg(validate_prior(set_prior("normal(0, 5)", class = "b"), f,
                             data = d))
  calls <- regmatches(m, gregexpr('set_prior[(].*?resp = "[^"]*"[)]', m,
                                   perl = TRUE))[[1]]
  expect_setequal(calls, c(
    'set_prior("normal(0, 5)", class = "b", resp = "y1")',
    'set_prior("normal(0, 5)", class = "b", dpar = "mub", resp = "cat")',
    'set_prior("normal(0, 5)", class = "b", dpar = "muc", resp = "cat")'))
  # and every one resolves
  for (cl in calls) {
    expect_length(ld_reached(f, d, eval(str2lang(cl))), 1L)
  }
})

test_that("a family brms does not have is refused by frmtmb's rule", {
  set.seed(42)
  n <- 200
  cl <- stats::rbinom(n, 1, 0.4)
  d <- data.frame(x = stats::rnorm(n))
  d$Y <- cbind(stats::rnorm(n, ifelse(cl == 1, 0, 3)),
               stats::rnorm(n, ifelse(cl == 1, 0, 4)))
  m <- ld_msg(validate_prior(set_prior("normal(0, 1)", class = "b"),
                             bf(Y ~ x) + mixture_mvn(K = 2, D = 2),
                             data = d))
  expect_match(m, "frmtmb's rule for a family", fixed = TRUE)
  expect_no_match(m, "as in brms", fixed = TRUE)
  expect_identical(ld_reached(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), d,
                              set_prior("normal(0, 1)", class = "b",
                                        dpar = "mu1d1")), "mu1d1_x")
})
