# Fix round 1, part B: the tests of F10 (us_t and diag_t are full rank)
# and F11 (an unstated full_rank leaves the check open) in
# tests/testthat/test-weights-unbounded.R. Rscript <this>
p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/tests/testthat/test-weights-unbounded.R"
txt <- paste(readLines(p), collapse = "\n")
i <- regexpr("\ntest_that(\"a registered structure is checked by its own rank\"",
             txt, fixed = TRUE)
stopifnot(i > 0)
txt <- substr(txt, 1L, i)
txt <- paste0(txt, '
test_that("a registered structure is checked by the rank it states", {
  # through 0.74.0 the check named the full-rank structures in a list,
  # so a structure another package registered was never one, and the
  # check failed open. It now reads the structure\'s full_rank. One that
  # states nothing leaves the check open: the covariance does not tell
  # whether a row\'s effect can reach any value (a latent with bounded
  # support keeps its maximum, below)
  local_frmtmb_registries()
  iid_lpdf <- function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
  }
  iid_vcov <- function(theta, blk) {
    matrix(exp(2 * theta[1]), 1, 1,
           dimnames = list(blk[["cnms"]], blk[["cnms"]]))
  }
  # an iid block of one coefficient under another name, full rank
  frmtmb_register_covstruct("myiid", lpdf = iid_lpdf, vcov = iid_vcov,
                            npar = 1L, sd_idx = 1L, full_rank = TRUE)
  # the same block, stating nothing
  frmtmb_register_covstruct("myiid2", lpdf = iid_lpdf, vcov = iid_vcov,
                            npar = 1L, sd_idx = 1L)
  # a block of two coefficients that moves only along (1, -1): its
  # covariance has rank one, and on columns a = b the effects add
  # nothing a row could use to reach its observation
  v <- c(1, -1) / sqrt(2)
  w <- c(1, 1) / sqrt(2)
  frmtmb_register_covstruct(
    "rankone",
    lpdf = function(b, theta, blk) {
      B <- RTMB::matrix(b, 2L)
      sum(RTMB::dnorm(as.vector(t(B) %*% v), 0, exp(theta[1]), log = TRUE)) +
        sum(RTMB::dnorm(as.vector(t(B) %*% w), 0, 1e-4, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(2 * theta[1]) * tcrossprod(v)
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = 1L, sd_idx = 1L)
  # a latent with bounded support, a scaled beta(2, 2) on (-0.3, 0.3):
  # positive definite at every parameter value, and no runaway, since a
  # row\'s effect cannot reach its observation
  frmtmb_register_covstruct(
    "boxre",
    lpdf = function(b, theta, blk) {
      z <- (b + 0.3) / 0.6
      sum(log(6 * z * (1 - z) / 0.6))
    },
    vcov = function(theta, blk) {
      matrix(0.3^2 / 5, 1, 1, dimnames = list(blk[["cnms"]], blk[["cnms"]]))
    },
    npar = 0L)
  d <- wu_data()
  d$w15 <- 1.5
  d$one <- 1
  d$one2 <- 1
  G <- Gamma(link = "log")
  # the stated full-rank block is refused as the built-in one is
  e_us <- tryCatch(frm(bf(y | weights(w15) ~ x + (1 | obs)), family = G,
                       data = d), error = conditionMessage)
  e_my <- tryCatch(frm(bf(y | weights(w15) ~ x + myiid(1 | obs)),
                       family = G, data = d), error = conditionMessage)
  expect_match(e_my, wu_refusal)
  expect_identical(sub("myiid(1 | obs)", "(1 | obs)", e_my, fixed = TRUE),
                   e_us)
  expect_true(frmtmb_covstruct_info("myiid")$full_rank)
  # the unstated one is not refused: the check is open for it, and the
  # reader says so
  expect_identical(frmtmb_covstruct_info("myiid2")$full_rank, NA)
  expect_false(frmtmb:::cs_full_rank(list(covstruct = "myiid2", dim = 1L,
                                          cnms = "Intercept",
                                          theta_idx = 1L)))
  # the bounded latent, unstated, fits and has its maximum
  rb <- wu_fit(bf(y | weights(w15) ~ x + boxre(1 | obs)), family = G,
               data = d)
  expect_true(is.finite(as.numeric(logLik(rb$fit))))
  expect_lt(exp(rb$fit$estimates$betad[1L]), 100)
  # the rank-one block, unstated, is not refused, and the model has its
  # maximum
  r <- wu_fit(bf(y | weights(w15) ~ x + rankone(0 + one + one2 | obs)),
              family = G, data = d)
  # the block adds nothing a row can use, so the maximum is the one
  # without it, to the 1e-4 sd its unused direction keeps; its own sd
  # is flat there, which the fit says in its warnings (wu_fit() keeps
  # them)
  r0 <- wu_fit(bf(y | weights(w15) ~ x), family = G, data = d)
  ll <- as.numeric(logLik(r$fit))
  expect_true(is.finite(ll))
  expect_equal(ll, as.numeric(logLik(r0$fit)), tolerance = 1e-6)
  expect_equal(unname(fixef(r$fit)[, 1]), unname(fixef(r0$fit)[, 1]),
               tolerance = 1e-4)
  expect_true(frmtmb_covstruct_info("us")$full_rank)
  expect_false(frmtmb_covstruct_info("rr")$full_rank)
})

test_that("a Student-t latent per row is full rank, and refused under
           weights that sum to more than the rows", {
  # through 0.74.0 us_t and diag_t were not full rank, so the check was
  # open for them: a weighted Gamma model with (1 | gr(obs, dist =
  # "student")) ran to a logLik of 4.7e7, the shape to 1e7. A t density
  # is positive everywhere, so a row\'s effect reaches its observation as
  # a gaussian one does
  d <- wu_data()
  d$w15 <- 1.5
  G <- Gamma(link = "log")
  expect_true(frmtmb_covstruct_info("us_t")$full_rank)
  expect_true(frmtmb_covstruct_info("diag_t")$full_rank)
  err <- tryCatch(frm(y | weights(w15) ~ x + (1 | gr(obs, dist = "student")),
                      family = G, data = d), error = conditionMessage)
  expect_match(err, wu_refusal)
  expect_match(err, "gr(obs, dist = \\"student\\")", fixed = TRUE)
  # the absent case: unweighted, the same model has its maximum
  r <- wu_fit(y ~ x + (1 | gr(obs, dist = "student")), family = G, data = d)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_lt(exp(r$fit$estimates$betad[1L]), 100)
})
')
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("ok\n")
