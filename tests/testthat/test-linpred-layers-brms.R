# brms 2.23.0 on the layers of test-linpred-layers.R, at the same parameter
# values (a Fixed_param brmsfit, brms_fixed_fit()). Pair p013's model and the
# others that a cs() term reaches. Seen to fail on 0.78.0 (rellib-r16):
# frm_curve() against posterior_linpred() read Inf, a length mismatch
# (dev/predread-cscurve.R).

skip_unless_brms()
skip_on_cran()
if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
  skip("set FRMTMB_BRMS_FIT_TESTS=true to run the fits against brms")
}

llb_data <- function() {
  set.seed(13)
  n <- 200
  d <- data.frame(x = round(stats::rnorm(n), 3), z = round(stats::rnorm(n), 3),
                  f = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- cut(d$x - 0.5 * d$z + stats::rlogis(n),
             c(-Inf, -0.8, 0.2, 1.2, Inf), labels = FALSE)
  d
}
llb_rel <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  expect_identical(length(a), length(b))
  max(abs(a - b) / pmax(1, abs(a), abs(b)))
}

test_that("frm_linpred, frm_curve and fitted match posterior_linpred", {
  d <- llb_data()
  fit <- frm(bf(y ~ x + cs(z)), family = sratio(), data = d)
  bfit <- brms_fixed_fit(brms::bf(y ~ x + cs(z)), brms::sratio(), d, fit)
  nd <- data.frame(x = c(-1, 0, 1), z = c(0.2, 0.2, -1))
  bl <- brms::posterior_linpred(bfit, newdata = nd, re_formula = NA)
  expect_identical(dim(bl)[2:3], c(3L, 3L))
  expect_lt(llb_rel(frm_linpred(fit, newdata = nd), bl[1, , ]), 1e-8)
  expect_lt(llb_rel(frm_curve(fit, newdata = nd,
                              simultaneous = FALSE)$.estimate, bl[1, , ]),
            1e-8)
  expect_lt(llb_rel(frm_lp_basis(fit, newdata = nd)$eta, bl[1, , ]), 1e-8)
  expect_lt(llb_rel(fitted(fit, newdata = nd, scale = "linear")[, "Estimate", ],
                    bl[1, , ]), 1e-8)
  # in sample
  expect_lt(llb_rel(frm_linpred(fit),
                    brms::posterior_linpred(bfit)[1, , ]), 1e-8)
  # dpar = "mu", transform = TRUE: the shared part, one number per row
  bm <- brms::posterior_linpred(bfit, newdata = nd, dpar = "mu")
  expect_lt(llb_rel(frm_linpred(fit, newdata = nd, dpar = "mu"), bm[1, ]),
            1e-8)
  bt <- brms::posterior_linpred(bfit, newdata = nd, transform = TRUE)
  expect_lt(llb_rel(frm_linpred(fit, newdata = nd, dpar = "mu",
                                type = "response"), bt[1, ]), 1e-8)
  bf <- fitted(bfit, newdata = nd, dpar = "mu", scale = "linear")
  expect_identical(dim(bf), c(3L, 4L))
  expect_lt(llb_rel(fitted(fit, newdata = nd, dpar = "mu",
                           scale = "linear")[, "Estimate"],
                    bf[, "Estimate"]), 1e-8)
})

test_that("a cs() factor and emmeans's rep.meas match", {
  skip_if_not_installed("emmeans")
  d <- llb_data()
  fit <- frm(bf(y ~ x + cs(f)), family = acat(), data = d)
  bfit <- brms_fixed_fit(brms::bf(y ~ x + cs(f)), brms::acat(), d, fit)
  nd <- data.frame(x = c(-1, 0, 1), f = factor(c("a", "b", "c"),
                                               levels = levels(d$f)))
  bl <- brms::posterior_linpred(bfit, newdata = nd, re_formula = NA)
  expect_lt(llb_rel(frm_curve(fit, newdata = nd,
                              simultaneous = FALSE)$.estimate, bl[1, , ]),
            1e-8)
  at <- list(x = c(-1, 1), f = "b")
  ea <- as.data.frame(emmeans::emmeans(fit, ~ x | rep.meas, at = at))
  eb <- suppressWarnings(as.data.frame(emmeans::emmeans(
    bfit, ~ x | rep.meas, at = at)))
  expect_identical(as.character(ea$rep.meas), as.character(eb$rep.meas))
  expect_lt(llb_rel(ea$emmean, eb$emmean), 1e-8)
  ea2 <- as.data.frame(emmeans::emmeans(fit, ~ x, at = at))
  eb2 <- suppressWarnings(as.data.frame(emmeans::emmeans(bfit, ~ x, at = at)))
  expect_lt(llb_rel(ea2$emmean, eb2$emmean), 1e-8)
})
