# refcat, the bare multinomial() and the simplex families in mvbf()
# against the Stan program brms 2.23.0 generates for them, at frmtmb's
# estimates: brms_lp_check() in helper-brms.R. With flat priors brms's
# log density is the log-likelihood, and its gradient vanishes at
# frmtmb's maximum. A refcat brms object is passed as brms spells it.
#
# A categorical category with its own formula is not here: brms 2.23.0
# stops on it in stancode() (brms-B in dev/famsimplex-findings.md),
# under refcat = "c" for the second predicted category too
# (dev/catmulti-findings.md).
#
# Stan compiles here, so the file is opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")

rcb_tight <- frmtmb_control(grad_tol = 1e-6, restarts = 3)

rcb_cat_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x, -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$f <- factor(apply(P, 1, function(p) {
    sample(c("a", "b", "c"), 1, prob = p)
  }), levels = c("a", "b", "c"))
  Yc <- t(vapply(seq_len(n), function(i) {
    as.numeric(stats::rmultinom(1, 10, P[i, ]))
  }, numeric(3)))
  colnames(Yc) <- c("y1", "y2", "y3")
  dd$Yc <- Yc
  dd$N <- 10L
  dd
}

test_that("categorical(refcat = 'b') is brms's likelihood", {
  skip_unless_brms_fit()
  dd <- rcb_cat_data(41)
  fit <- frm(bf(f ~ x), family = brms::categorical(refcat = "b"),
             data = dd, control = rcb_tight)
  expect_identical(names(fixef_by_dpar(fit)), c("mua", "muc"))
  brms_lp_check(brms::bf(f ~ x), brms::categorical(refcat = "b"), dd, fit)
})

test_that("multinomial() and multinomial(refcat = 'y2') are brms's", {
  skip_unless_brms_fit()
  dd <- rcb_cat_data(43, n = 200)
  fit <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dd,
             control = rcb_tight)
  brms_lp_check(brms::bf(Yc | trials(N) ~ x), brms::multinomial(), dd, fit)
  f2 <- frm(bf(Yc | trials(N) ~ x), family = multinomial(refcat = "y2"),
            data = dd, control = rcb_tight)
  brms_lp_check(brms::bf(Yc | trials(N) ~ x),
                brms::multinomial(refcat = "y2"), dd, f2)
})

test_that("refcat = NA is brms's likelihood, along its flat direction", {
  skip_unless_brms_fit()
  dd <- rcb_cat_data(44, n = 200)
  # the maximum is a line, and brms's gradient vanishes anywhere on it
  allow_warnings(
    fit <- frm(bf(Yc | trials(N) ~ x), family = multinomial(refcat = NA),
               data = dd, control = rcb_tight),
    c("Standard errors are not available",
      "Large maximum absolute gradient"))
  brms_lp_check(brms::bf(Yc | trials(N) ~ x),
                brms::multinomial(refcat = NA), dd, fit)
})

test_that("dirichlet() inside mvbf() is brms's likelihood", {
  skip_unless_brms_fit()
  set.seed(45)
  n <- 300
  dd <- data.frame(x = stats::rnorm(n))
  E <- cbind(0, 0.4 + 0.5 * dd$x, -0.3 - 0.4 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  G <- matrix(stats::rgamma(n * 3, shape = 6 * P), n, 3)
  Y <- G / rowSums(G)
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
  fit <- frm(mvbf(bf(Y ~ x, family = dirichlet()),
                  bf(w ~ x, family = gaussian())), data = dd,
             control = rcb_tight)
  brms_lp_check(brms::bf(Y ~ x, family = brms::dirichlet()) +
                  brms::bf(w ~ x, family = gaussian()) +
                  brms::set_rescor(FALSE),
                list(brms::dirichlet(), gaussian()), dd, fit)
})

test_that("logistic_normal() inside mvbf() is brms's model block", {
  skip_unless_brms_fit()
  set.seed(46)
  n <- 300
  dd <- data.frame(x = stats::rnorm(n))
  W <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .4, .4, 1), 2))
  E <- cbind(0, 0.4 + 0.5 * dd$x + 0.6 * W[, 1],
             -0.3 - 0.4 * dd$x + 0.9 * W[, 2])
  Y <- exp(E) / rowSums(exp(E))
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
  fit <- frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                  bf(w ~ x, family = gaussian())), data = dd,
             control = rcb_tight)
  # brms 2.23.0 cannot compile this model: its generated quantities read
  # Llncor and nlncor where it declared Llncor_Y and nlncor_Y
  # (dev/catmulti-out/brms-fits.txt). The repair touches the generated
  # quantities alone, which log_prob() does not evaluate.
  fix <- function(code) {
    out <- gsub("multiply_lower_tri_self_transpose(Llncor);",
                "multiply_lower_tri_self_transpose(Llncor_Y);", code,
                fixed = TRUE)
    out <- gsub("[nlncor] lncor;", "[nlncor_Y] lncor;", out, fixed = TRUE)
    expect_false(identical(out, as.character(code)))
    out
  }
  brms_lp_check(brms::bf(Y ~ x, family = brms::logistic_normal()) +
                  brms::bf(w ~ x, family = gaussian()) +
                  brms::set_rescor(FALSE),
                list(brms::logistic_normal(), gaussian()), dd, fit,
                code_fix = fix)
})

test_that("a Y_a response maps onto brms's Llncor_Ya", {
  skip_unless_brms_fit()
  set.seed(47)
  n <- 300
  dd <- data.frame(x = stats::rnorm(n))
  W <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .4, .4, 1), 2))
  E <- cbind(0, 0.4 + 0.5 * dd$x + 0.6 * W[, 1],
             -0.3 - 0.4 * dd$x + 0.9 * W[, 2])
  Y <- exp(E) / rowSums(exp(E))
  colnames(Y) <- c("a", "b", "c")
  dd$Y_a <- Y
  dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
  fit <- frm(mvbf(bf(Y_a ~ x, family = logistic_normal()),
                  bf(w ~ x, family = gaussian())), data = dd,
             control = rcb_tight)
  # brms-C again, under brms's Stan name of the response
  fix <- function(code) {
    out <- gsub("multiply_lower_tri_self_transpose(Llncor);",
                "multiply_lower_tri_self_transpose(Llncor_Ya);", code,
                fixed = TRUE)
    out <- gsub("[nlncor] lncor;", "[nlncor_Ya] lncor;", out, fixed = TRUE)
    expect_false(identical(out, as.character(code)))
    out
  }
  brms_lp_check(brms::bf(Y_a ~ x, family = brms::logistic_normal()) +
                  brms::bf(w ~ x, family = gaussian()) +
                  brms::set_rescor(FALSE),
                list(brms::logistic_normal(), gaussian()), dd, fit,
                code_fix = fix)
})
