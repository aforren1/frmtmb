# A 2 x 2 correlation matrix is one chart coordinate, so confint(parm =)
# takes its correlation's name, as it takes sd_<group>__<term>; 0.71.0
# refused "lncor__b__c" and "rescor__y1__y2" with "does not stand for a
# single internal one", and advised hypothesis(method = 'profile'),
# which then failed ("not linear"). Also drop1() on a family without
# mu, which 0.71.0 refused with "No linear predictor for dpar 'mu'".

c2_ln_data <- function(seed, n = 300, K = 3) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  D <- K - 1L
  R <- matrix(0.4, D, D)
  diag(R) <- 1
  W <- matrix(stats::rnorm(D * n), n) %*% chol(R)
  E <- cbind(0, sapply(seq_len(D), function(k) 0.3 * k + 0.5 * dd$x) + W)
  Y <- exp(E) / rowSums(exp(E))
  colnames(Y) <- letters[seq_len(K)]
  dd$Y <- Y
  dd
}

test_that("confint(parm =) takes a 2 x 2 correlation's name", {
  dd <- c2_ln_data(1)
  fit <- frm(bf(Y ~ x), family = logistic_normal(), data = dd)
  expect_message(ci <- confint(fit, parm = "lncor__b__c"),
                 "is lncor_raw_1")
  expect_identical(ci, confint(fit, parm = "lncor_raw_1"))
  # rescor of two responses
  set.seed(2)
  d2 <- data.frame(x = stats::rnorm(200))
  e <- matrix(stats::rnorm(400), 200) %*% chol(matrix(c(1, .5, .5, 1), 2))
  d2$y1 <- 1 + d2$x + e[, 1]
  d2$y2 <- -1 + 0.5 * d2$x + e[, 2]
  fr <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE), data = d2)
  expect_message(cr <- confint(fr, parm = "rescor__y1__y2"),
                 "is thetar_1")
  expect_identical(cr, confint(fr, parm = "thetar_1"))
})

test_that("a wider correlation block is refused with advice that works", {
  dd <- c2_ln_data(3, K = 4)
  fit <- frm(bf(Y ~ x), family = logistic_normal(), data = dd)
  expect_error(confint(fit, parm = "lncor__b__c"),
               "hypothesis(fit, \"lncor__b__c = 0\", class = NULL)",
               fixed = TRUE)
  h <- hypothesis(fit, "lncor__b__c = 0", class = NULL)
  expect_true(is.finite(h$hypothesis$Estimate))
  expect_true(is.finite(h$hypothesis$Est.Error))
})

test_that("drop1() drops a term from every category's predictor", {
  set.seed(4)
  n <- 300
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x, -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$f <- factor(apply(P, 1, function(p) sample(letters[1:3], 1, prob = p)))
  fit <- frm(bf(f ~ x + z), family = categorical(), data = dd)
  tab <- drop1(fit, test = "Chisq")
  expect_identical(rownames(tab), c("<none>", "x", "z"))
  # each term leaves both logits: 2 degrees of freedom, and the LRT is
  # the refit's
  expect_identical(tab$Df[-1L], c(2, 2))
  f_x <- frm(bf(f ~ z), family = categorical(), data = dd)
  lrt <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(f_x)))
  expect_lt(abs(tab["x", "LRT"] - lrt), 1e-8 * abs(lrt))
  # every location predictor with its own formula: nothing to drop from
  fo <- frm(bf(f ~ x, mub ~ z, muc ~ z), family = categorical(), data = dd)
  expect_error(drop1(fo), "has a formula of its own")
})
