# The silent wrong answers the ported brms suite found, and the refusal
# or the answer each one now gets. Lane wt-adefects; the constructions
# and the before-and-after numbers are in dev/adefects-findings.md, and
# dev/adefects-repro-base.R reproduces every one of them on the base
# build (frmtmb 0.60.0) beside brms 2.23.0.
#
# Every block below was SEEN FAILING against that base build
# (dev/adefects-log/newtests-base.txt). Where the defect was a wrong
# NUMBER rather than a missing refusal, the assertion is on the number:
# a block that only asserted `expect_error` would pass on a build that
# raises the right error for the wrong reason.

ad_data <- function(seed = 20260917) {
  set.seed(seed)
  d <- data.frame(g1 = rep(c(1, 1, 2, 2), each = 6),
                  g2 = rep(c(1, 2, 1, 2), each = 6),
                  t = c(1:6, 1:6, 1:6, 7:12))
  d$x <- rnorm(24)
  d$g <- interaction(d$g1, d$g2)
  d$gf1 <- factor(d$g1)
  d$gf2 <- factor(d$g2)
  d$y <- d$x + as.numeric(stats::filter(rnorm(24), 0.6, "recursive"))
  d
}

ad_term <- function(txt) {
  stats::as.formula(paste("y ~ x +", txt))
}

# ---------------------------------------------------------------- D1

test_that("an autocorrelation term refuses an expression as its time index", {
  d <- ad_data()
  # the control still fits; its value is the reference every comparison
  # below is made against, rather than a number typed into the test
  ctl <- frm(y ~ x + ar(t, g, cov = TRUE), d)
  expect_true(is.finite(as.numeric(logLik(ctl))))

  # the defect: the base build fitted `x + t` AS the time index and
  # reported a different likelihood with nothing said
  for (tm in c("ar(x + t, g, cov = TRUE)", "ar(t - 10 * x, g, cov = TRUE)",
               "ma(x + t, g, cov = TRUE)", "arma(x + t, g, cov = TRUE)",
               "cosy(x + t, g)", "unstr(x + t, g)")) {
    expect_error(frm(ad_term(tm), d), class = "frmtmb_error", label = tm)
    expect_error(frm(ad_term(tm), d), "ONE variable name", label = tm)
  }
  # a literal and a string are not variable names either, as in brms
  expect_error(frm(y ~ x + ar(1, g, cov = TRUE), d), "ONE variable name")
  expect_error(frm(y ~ x + ar("t", g, cov = TRUE), d), "ONE variable name")
  # the spellings brms accepts still fit, and to the control's value
  expect_equal(as.numeric(logLik(frm(y ~ x + ar(t, gr = g, cov = TRUE), d))),
               as.numeric(logLik(ctl)), tolerance = 1e-8)
  expect_s3_class(frm(y ~ x + ar(gr = g, cov = TRUE), d), "frmtmb_fit")
  expect_s3_class(frm(y ~ x + ar(NA, g, cov = TRUE), d), "frmtmb_fit")
})

# ---------------------------------------------------------------- D2

test_that("gr = takes variable names crossed by `:` and nothing else", {
  d <- ad_data()
  ctl <- frm(y ~ x + ar(t, g, cov = TRUE), d)

  # `gr = a:b` CROSSES, on numeric codes as well as on factors. R's `:`
  # is the sequence operator for numbers, so the base build read
  # `gr = g1:g2` as a length-one vector and reported repeated times;
  # the two spellings now agree with each other and with the
  # interaction column the control groups by.
  num <- frm(y ~ x + ar(t, gr = g1:g2, cov = TRUE), d)
  fac <- frm(y ~ x + ar(t, gr = gf1:gf2, cov = TRUE), d)
  expect_equal(as.numeric(logLik(num)), as.numeric(logLik(ctl)),
               tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fac)), as.numeric(logLik(ctl)),
               tolerance = 1e-8)
  # three crossed names are legal grammar, as they are in brms
  expect_s3_class(frm(y ~ x + ar(t, gr = g1:g2:g, cov = TRUE), d),
                  "frmtmb_fit")

  # every other operator is refused rather than evaluated. `g1/g2` was
  # the silent one: series (1, 1) and (2, 2) both have ratio 1 and
  # merged, for a likelihood that differed in the fifth decimal.
  for (gr in c("g1/g2", "g1 + g2", "g1 * g2", "interaction(g1, g2)",
               "factor(g)")) {
    tm <- paste0("ar(t, gr = ", gr, ", cov = TRUE)")
    expect_error(frm(ad_term(tm), d), class = "frmtmb_error", label = gr)
    expect_error(frm(ad_term(tm), d), "combined by", label = gr)
  }
})

# ---------------------------------------------------------------- D3

test_that("a hypothesis with no relation is refused, as in brms", {
  set.seed(20260917)
  dh <- data.frame(y = rnorm(40), Age = rnorm(40), Trt = rnorm(40))
  fh <- frm(y ~ Age + Trt, dh)

  # the defect: `"Age"` was answered as the row of `"Age = 0"`
  expect_error(hypothesis(fh, "Age"), class = "frmtmb_error")
  expect_error(hypothesis(fh, "Age"), "states no relation")
  expect_error(hypothesis(fh, "Age + Trt"), "states no relation")
  # one bad string in a vector refuses the call
  expect_error(hypothesis(fh, c("Age = 0", "Trt")), "states no relation")

  # the relation spellings still answer, and with the same numbers
  h0 <- hypothesis(fh, "Age = 0")$hypothesis
  expect_equal(h0$Estimate, unname(fixef(fh)$mu[["Age"]]), tolerance = 1e-10)
  expect_true(is.finite(hypothesis(fh, "Age > 0")$hypothesis$Est.Error))
  expect_true(is.finite(hypothesis(fh, "Age < Trt")$hypothesis$Est.Error))
  expect_true(is.finite(hypothesis(fh, "exp(Age) = 1")$hypothesis$Est.Error))
})

# ---------------------------------------------------------------- D4

test_that("fit$data is the model frame, not a partial match of data2", {
  set.seed(20260917)
  A <- diag(6)
  A[A == 0] <- 0.3
  dimnames(A) <- list(1:6, 1:6)
  dd <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
  dd$y <- dd$x + rnorm(6)[dd$g] + rnorm(30)
  fa <- frm(y ~ x + (1 | gr(g, cov = A)), dd, data2 = list(A = A))

  # the defect: with no `data` element, `$` partial-matched `data2`, so
  # a script written for brms read the covariance matrix rather than the
  # 30-row frame, and nothing said so
  expect_true("data" %in% names(fa))
  expect_s3_class(fa$data, "data.frame")
  expect_identical(fa[["data", exact = TRUE]], fa$data)
  expect_false(identical(fa$data, fa$data2))
  expect_identical(names(fa$data2), "A")

  # it is the model frame, which is what model.frame() already
  # returned. NOT brms's `data`, which is the validated RAW data: on a
  # model with a transformed term the two differ, and the block below
  # pins that difference rather than claiming a parity the package does
  # not have.
  expect_identical(fa$data, model.frame(fa))
  expect_equal(nrow(fa$data), 30L)
  expect_true(all(c("y", "x", "g") %in% names(fa$data)))

  # A REAL element, not a field served by a `$` method. These three
  # fail on a design that keeps one copy and serves `data` from the
  # frame through `$.frmtmb_fit`, which is the alternative punch round 1
  # weighed: `names()`, `[[` and `str()` would not show it there.
  expect_true("data" %in% names(fa))
  expect_false(is.null(fa[["data", exact = TRUE]]))
  expect_null(getS3method("$", "frmtmb_fit", optional = TRUE))

  # SHARED in memory, and the shared-ness is what value equality cannot
  # see. The control below is a COPY of the same frame: `data.frame()`
  # drops the `terms` attribute, so it is not `identical()` to the
  # original, but it is a distinct object at a distinct address, which
  # is what the address comparison has to tell apart.
  addr <- function(x) {
    strsplit(trimws(utils::capture.output(.Internal(inspect(x)))[1L]),
             "[[:space:]]+")[[1L]][1L]
  }
  expect_identical(addr(fa$data), addr(fa$frame[["data_frame"]]))
  expect_false(identical(addr(fa$data), addr(data.frame(fa$data))))

  # and NOT shared ON DISK, because R's serializer does not deduplicate
  # a shared value. The lane first claimed the element was free; it
  # costs a frame in every saved fit, 5.3% of the raw bytes and 14.0%
  # of the gzipped ones on a 20,000-row fit
  # (dev/adefects-log/p1-dollar.txt). Under a design that keeps one copy
  # and serves `data` through `$.frmtmb_fit`, this delta is 0.
  #
  # Measured as a RATIO to the COLUMNS this run serializes, never as a
  # byte count, and to the columns rather than to the frame: a model
  # frame carries a `terms` attribute whose environment is the caller's,
  # and serializing that drags the whole calling frame in, which read
  # 3,696,708 bytes here against 911 for the columns alone.
  no_data <- fa
  no_data$data <- NULL
  cols <- length(serialize(lapply(fa$data, identity), NULL))
  delta <- length(serialize(fa, NULL)) - length(serialize(no_data, NULL))
  expect_gt(delta, 0.9 * cols)
  expect_lt(delta, 3 * cols)

  # a model with a transformed term: the frame carries the TERM, where
  # brms's `data` carries the raw column and every unused one
  dd$z <- rnorm(30)
  fo <- frm(y ~ x + offset(z), dd)
  expect_true("offset(z)" %in% names(fo$data))
  expect_false("z" %in% names(fo$data))

  # the other brmsfit field names read NULL rather than a neighbour:
  # `data` was the one that collided. `exact = FALSE` is `$`'s own
  # partial matching, which is the read being tested; plain `[[` is
  # exact and would pass on a build that still has the defect.
  for (nm in c("formula", "family", "ranef", "criteria", "version",
               "algorithm", "backend", "basis", "stanvars", "model")) {
    expect_null(fa[[nm, exact = FALSE]], label = nm)
  }
  # the inverse, so that NULL above is not NULL of everything: the two
  # names that DO partial-match still do
  expect_identical(fa[["estimat", exact = FALSE]], fa$estimates)
  expect_identical(fa[["data2", exact = FALSE]], fa$data2)
})

# ---------------------------------------------------------------- D6

test_that("newdata may omit a grouping column when new levels are allowed", {
  set.seed(20260917)
  dg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
  dg$y <- dg$x + rnorm(8)[dg$g] + rnorm(40)
  fg <- frm(y ~ x + (1 | g), dg)

  nd <- data.frame(x = c(0, 1))
  nd_new <- data.frame(x = c(0, 1), g = factor(c("99", "98")))

  # the defect: this stopped at base R's "object 'g' not found", where
  # brms's validate_newdata() fills the column with NA and answers
  # (measured on brmsfit_example1 with `visit` removed)
  p_fill <- predict(fg, newdata = nd, allow_new_levels = TRUE)
  p_new <- predict(fg, newdata = nd_new, allow_new_levels = TRUE)
  p_pop <- predict(fg, newdata = nd, re_formula = NA)
  # an absent column means an unseen level, and an unseen level is the
  # population value for a maximum-likelihood fit
  expect_equal(unname(p_fill), unname(p_new), tolerance = 1e-12)
  expect_equal(unname(p_fill), unname(p_pop), tolerance = 1e-12)

  # and an unseen level carries the block's variance, which the
  # population prediction does not
  s_fill <- predict(fg, newdata = nd, allow_new_levels = TRUE,
                    se.fit = TRUE)$se.fit
  s_new <- predict(fg, newdata = nd_new, allow_new_levels = TRUE,
                   se.fit = TRUE)$se.fit
  s_pop <- predict(fg, newdata = nd, re_formula = NA, se.fit = TRUE)$se.fit
  expect_equal(unname(s_fill), unname(s_new), tolerance = 1e-12)
  expect_true(all(s_fill > s_pop))

  # without allow_new_levels it is refused BY NAME. brms stops here too,
  # on base R's own error; this one is classed and says what to pass.
  expect_error(predict(fg, newdata = nd), class = "frmtmb_error")
  expect_error(predict(fg, newdata = nd), "allow_new_levels = TRUE")
  expect_error(predict(fg, newdata = nd), "`g`", fixed = TRUE)
  # a column that IS there keeps the refusal it had: a level the fit
  # never saw is still named by the old message
  expect_error(predict(fg, newdata = nd_new), "New levels in grouping")
})

# ---------------------------------------------------------------- D7

test_that("log_lik() on a maximum-likelihood fit says to sample", {
  set.seed(20260917)
  dd <- data.frame(x = rnorm(40))
  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
  fit <- frm(y ~ x, dd)

  # the defect: frmtmb defined no log_lik generic, so this was base R's
  # "could not find function", where loo() and waic() name the route
  expect_true(is.function(frmtmb::log_lik))
  expect_error(log_lik(fit), class = "frmtmb_error")
  expect_error(log_lik(fit), "Sample first")
  # the same shape as the two refusals it is modelled on
  expect_error(loo(fit), class = "frmtmb_error")
  expect_error(waic(fit), class = "frmtmb_error")

  # registered on the OWNER's generic as well as on frmtmb's, so the
  # refusal is still reached once rstantools is loaded and its generic
  # takes the binding
  ns <- parseNamespaceFile("frmtmb", dirname(find.package("frmtmb")))
  m <- ns$S3methods
  row <- m[, 1] == "log_lik" & m[, 2] == "frmtmb_fit"
  expect_equal(sum(row), 2L)
  expect_true(any(is.na(m[row, 4])))
  expect_true(any(m[row, 4] %in% "rstantools"))
  # and the owner table knows who owns the name
  own <- get("frm_generic_owners", envir = asNamespace("frmtmb"))
  expect_identical(own[["log_lik"]], "rstantools")
})
