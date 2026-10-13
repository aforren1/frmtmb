# A data column or a nonlinear parameter NAMED like a special is a plain
# variable (punch round 2, r3). The parser found the specials with
# all.names(), which returns symbols as well as calls, so a column named
# mo, mi, ar, ma, arma, cosy or unstr was read as the special: on 0.69.0
# and on the lane before the fix, y ~ mo * z was refused as a crossed
# mo() term, y ~ ar * z as a crossed autocorrelation term, and a
# nonlinear parameter named ma as a residual correlation term in the
# body (dev/moparity-rev2-moname.R, dev/moparity-r3-names.R). Seen to
# fail on the lane build of punch round 1.

special_named_data <- function() {
  set.seed(2)
  n <- 200
  d <- data.frame(v = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(sample(1:15, n, TRUE)))
  d$y <- 1 + 0.3 * d$v * d$z + stats::rnorm(15)[d$g] + stats::rnorm(n)
  d$yp <- exp(0.5 + 0.4 * d$z) + stats::rnorm(n, 0, 0.2)
  d
}

same_fit <- function(fo, d, nm) {
  ref <- allow_warnings(frm(bf(fo), family = gaussian(), data = d),
                        c("Boundary", "singular"))
  names(d)[names(d) == "v"] <- nm
  fo2 <- stats::as.formula(gsub("\\bv\\b", nm, deparse1(fo)))
  got <- allow_warnings(frm(bf(fo2), family = gaussian(), data = d),
                        c("Boundary", "singular"))
  ll <- as.numeric(logLik(ref))
  expect_lt(abs(as.numeric(logLik(got)) - ll), 1e-10 * abs(ll),
            label = paste(nm, deparse1(fo2)))
}

test_that("a column named like a special is a plain variable", {
  d <- special_named_data()
  for (nm in c("mo", "mi", "ar", "ma", "arma", "cosy", "unstr")) {
    same_fit(y ~ v * z, d, nm)
    same_fit(y ~ v:z, d, nm)
    same_fit(y ~ v + (1 + v || g), d, nm)
    same_fit(y ~ v + (1 + v | g), d, nm)
  }
})

test_that("a nonlinear parameter named like a special is a parameter", {
  d <- special_named_data()
  pr <- set_prior("normal(0, 5)", nlpar = "b")
  ref <- frm(bf(yp ~ a + exp(b * z), a ~ 1, b ~ 1, nl = TRUE),
             family = gaussian(), data = d, prior = pr)
  ll <- as.numeric(logLik(ref))
  for (nm in c("ar", "ma", "arma", "cosy", "unstr")) {
    fo <- bf(stats::as.formula(sprintf("yp ~ %s + exp(b * z)", nm)),
             stats::as.formula(sprintf("%s ~ 1", nm)), b ~ 1, nl = TRUE)
    got <- frm(fo, family = gaussian(), data = d, prior = pr)
    expect_lt(abs(as.numeric(logLik(got)) - ll), 1e-10 * abs(ll),
              label = nm)
  }
})

test_that("the specials themselves are still read as specials", {
  # the guards absent from the cases above still fire on real calls
  d <- special_named_data()
  d$x <- sample(0:3, nrow(d), TRUE)
  d$w <- sample(0:2, nrow(d), TRUE)
  expect_error(frm(bf(y ~ mo(x):mo(w)), family = gaussian(), data = d),
               "cannot interact with another")
  expect_error(frm(bf(y ~ z * ar(gr = g)), family = gaussian(), data = d),
               "autocorrelation term cannot be crossed")
  expect_error(frm(bf(yp ~ a * ar(z), a ~ 1, nl = TRUE), family = gaussian(),
                   data = d),
               "Residual correlation terms are not supported")
  fit <- frm(bf(y ~ mo(x) * z), family = gaussian(), data = d)
  expect_true("mox:z" %in% rownames(fixef(fit)))
})
