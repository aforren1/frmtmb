test_that("a bound typed at a computed end of the interval is that end", {
  # the 0.76.0 release review, finding 3: a ring of 20 with chords,
  # row-standardized, has its largest eigenvalue at 1, and 1 / it
  # computes to 0.999999999999998 or to 1 by how LAPACK rounds; ub = 1
  # was refused on the first, as "outside ... (-1.72756, 1)"
  n <- 20
  A <- matrix(0, n, n)
  for (i in 1:n) {
    j <- i %% n + 1
    A[i, j] <- A[j, i] <- 1
    j <- (i + 1) %% n + 1
    A[i, j] <- A[j, i] <- 1
  }
  A[1, 10] <- A[10, 1] <- 1
  W <- A / rowSums(A)
  set.seed(3)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- as.vector(solve(diag(n) - 0.4 * W, 1 + 0.5 * d$x +
                           stats::rnorm(n, 0, 0.7)))
  d2 <- list(W = W)
  for (type in c("lag", "error")) {
    fo <- stats::as.formula(paste0("y ~ x + sar(W, type = '", type, "')"))
    cl <- paste0(type, "sar")
    f0 <- frm(fo, data = d, data2 = d2)
    dp <- default_prior(fo, data = d, data2 = d2)
    hi <- dp$ub[dp$class == cl]
    # the end the table lists is the computed one: 1 within its ulps
    expect_lt(abs(hi - 1), 64 * .Machine$double.eps)
    p1 <- set_prior("", class = cl, ub = 1)
    vp <- validate_prior(p1, fo, data = d, data2 = d2)
    expect_true(cl %in% vp$class)
    f1 <- frm(fo, data = d, data2 = d2, prior = p1)
    expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
    # the guard absent: past the end by more than the ulps an end
    # carries, the bound is refused, and the message prints the end
    e <- tryCatch(frm(fo, data = d, data2 = d2,
                      prior = set_prior("", class = cl, ub = 1 + 1e-9)),
                  error = identity)
    expect_s3_class(e, "frmtmb_error")
    expect_match(conditionMessage(e), "is outside the interval of sar(W",
                 fixed = TRUE)
    expect_match(conditionMessage(e), format(hi, digits = 15), fixed = TRUE)
  }
})
