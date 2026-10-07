# fitted() on an ordinal fit whose disc predictor has no fixed column
# (lane surface, 2026-10-06; found by the gpby review).
#
# Seen to fail on 0.68.1 (rellib-r6): fitted() stopped with "requires
# numeric/complex matrix/vector arguments" in lp_eta_design(), in
# sample and on newdata, because the estimates carry no betad when no
# disc column is estimated (dev/surface-repros.R item 8,
# dev/gpby-p1-disc2.R).

disc_data <- function() {
  set.seed(11)
  n <- 120
  d <- data.frame(x = round(stats::runif(n, 0, 5), 1), z = stats::rnorm(n),
                  g = factor(sample(letters[1:8], n, TRUE)))
  lat <- sin(d$x) + 0.5 * d$z + stats::rlogis(n)
  d$y <- factor(cut(lat, stats::quantile(lat, c(0, 0.3, 0.6, 1)),
                    include.lowest = TRUE, labels = FALSE), ordered = TRUE)
  d
}

for (rhs in c("0 + (1 | g)", "0 + gp(x, k = 6)")) {
  test_that(paste("fitted() runs on disc ~", rhs), {
    d <- disc_data()
    fit <- suppressWarnings(frm(bf(y ~ z, stats::as.formula(paste("disc ~",
                                                               rhs))),
                                data = d, family = cumulative()))
    expect_null(fit$estimates[["betad"]])
    p <- allow_warnings(fitted(fit), "could not be recovered")
    expect_equal(dim(p), c(nrow(d), 4L, 3L))
    # each row's category probabilities sum to one
    expect_equal(unname(rowSums(p[, "Estimate", ])), rep(1, nrow(d)))
    pn <- allow_warnings(fitted(fit, newdata = d[1:5, ]),
                         "could not be recovered")
    expect_equal(pn[, "Estimate", ], p[1:5, "Estimate", ])
    # the fixed part of a predictor with no fixed column is zero, as the
    # objective forms it: the disc predictor is the group or gp() part
    lp <- fit$frame$linpreds[["disc"]] %||%
      fit$frame$linpreds[[grep("disc", names(fit$frame$linpreds))[1L]]]
    expect_identical(frmtmb:::lp_fixed_part(lp$X, fit$estimates, lp),
                     numeric(nrow(d)))
  })
}

test_that("a predictor with a fixed column is unchanged (control)", {
  d <- disc_data()
  fit <- suppressWarnings(frm(bf(y ~ z, disc ~ 1 + (1 | g)), data = d,
                              family = cumulative()))
  lp <- fit$frame$linpreds[[grep("disc", names(fit$frame$linpreds))[1L]]]
  expect_identical(frmtmb:::lp_fixed_part(lp$X, fit$estimates, lp),
                   drop(as.matrix(lp$X %*% fit$estimates[["betad"]][
                     lp$idx])))
  expect_equal(dim(fitted(fit)), c(nrow(d), 4L, 3L))
})
