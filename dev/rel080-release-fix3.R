# The release review of 0.80.0 (dev/reviews/2026-11-07-release.md, verdict
# FIX), applied as its items word them, after dev/rel080-release-fix2.R:
#
#   Rscript dev/rel080-release-fix3.R
#
# NEWS (items 2, 3, 4) is written by dev/rel080-news.R and the record
# (items 5, 8) by hand; this script edits code, tests, the checklist, the
# survey's generator, the commit message and the backlog. Each old text
# must occur exactly once.
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
rd <- function(p) paste(readLines(file.path(REL, p), warn = FALSE),
                        collapse = "\n")
wr <- function(x, p) {
  con <- file(file.path(REL, p), open = "wb")
  on.exit(close(con))
  writeLines(strsplit(x, "\n", fixed = TRUE)[[1L]], con, sep = "\r\n")
}
edit <- function(p, old, new, what) {
  x <- rd(p)
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(p, ": ", what, ": the old text occurs ", n, " times")
  x <- sub(old, new, x, fixed = TRUE)
  if (grepl("[.]R$", p)) invisible(parse(text = x, keep.source = FALSE))
  wr(x, p)
  cat("EDIT", p, "|", what, "\n")
}
append_to <- function(p, add, what) {
  x <- rd(p)
  x <- paste0(sub("\n*$", "", x), "\n\n", add)
  if (grepl("[.]R$", p)) invisible(parse(text = x, keep.source = FALSE))
  wr(x, p)
  cat("APPEND", p, "|", what, "\n")
}

# item 1 (MAJOR): no body-basis descent on the refit of a probe's winner
edit("R/fit.R",
"    if (is.null(bw)) {
      opt <- nl_basis_descent(obj, opt, frame, ctl_opt, bounds, par_units,
                              tally, vb)
    }",
"    # not on the refit of a sign probe's winner: from there a second
    # descent in the body's Jacobian basis ran on along a ridge to the
    # linear limit, to a point whose standard errors read finite, and the
    # fit lost the warning that it is not a maximum (binomial
    # Michaelis-Menten, 4.6 short; release review of 0.80.0, item 1)
    if (is.null(bw) && !isTRUE(attr(start, \"frmtmb_nl_probe\"))) {
      opt <- nl_basis_descent(obj, opt, frame, ctl_opt, bounds, par_units,
                              tally, vb)
    }",
"item 1, fit_assembled()")
append_to("tests/testthat/test-nl-start.R",
"test_that(\"a probe winner on a ridge is not made silent by a second descent\", {
  # binomial Michaelis-Menten (review of lane nlstart, fresh row mm 22):
  # the sign probe's winner runs to K -> Inf along the ridge to the linear
  # limit, 4.6 below the maximum, and its standard errors say so. A second
  # descent in the body's Jacobian basis from there slid further along
  # the ridge, the standard errors read finite and the fit was silent
  set.seed(22001)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  # the draws of the review's design, in its order
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  d$nt <- 10L
  d$y <- rbinom(n, 10, plogis(3 * d$xp / (1 + d$xp) - 2))
  fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm ~ 1, K ~ 1, nl = TRUE)
  got <- caught(frm(fm, family = binomial(), data = d))
  at <- caught(frm(fm, family = binomial(), data = d, se = FALSE,
                   start = list(beta = c(\"Vm_(Intercept)\" = -0.1,
                                         \"K_(Intercept)\" = -6))))
  l_def <- as.numeric(logLik(got$value))
  l_max <- as.numeric(logLik(at$value))
  warned <- any(vapply(got$conds, inherits, NA, \"warning\"))
  # the default fit reaches the maximum, or it says that it may not
  expect_true(l_max - l_def <= 1e-6 * abs(l_max) || warned)
})",
"item 1, the reviewer's test")

# item 2 (MAJOR): `0 + f` designs centered outside the span of the constant
edit("R/frame.R",
"#' A fixed design with every column but the intercept centered, as
#' opt_center() centers it (the design itself where there is no
#' intercept). A column whose centered root mean square is at most",
"#' A fixed design with every column but the intercept centered, as
#' opt_center() centers it (the design itself where no column is the
#' intercept and the columns do not span the constant; where they span
#' it, as `0 + f` does, every column outside that span is centered). A
#' column whose centered root mean square is at most",
"item 2, design_center()'s roxygen")
edit("R/frame.R",
"design_center <- function(X) {
  ic <- match(\"(Intercept)\", colnames(X))
  if (is.na(ic)) return(X)
  m <- colMeans(X)
  m[ic] <- 0
  Xc <- X - rep(m, each = nrow(X))
  raw <- sqrt(colSums(X^2))
  cen <- sqrt(colSums(Xc^2))
  flat <- seq_len(ncol(X)) != ic & !(cen > design_alias_tol * raw)
  Xc[, flat] <- 0
  Xc
}",
"design_center <- function(X) {
  ic <- match(\"(Intercept)\", colnames(X))
  stay <- rep(FALSE, ncol(X))
  if (is.na(ic)) {
    # no intercept column, but columns that span the constant (`0 + f`):
    # a multiple of the constant taken from a column outside that span
    # leaves the column space as it was, and a column far from zero is as
    # collinear with the constant there as beside an intercept (release
    # review of 0.80.0, item 2)
    one <- rep(1, nrow(X))
    q <- qr(X, tol = design_alias_tol)
    if (!(sqrt(sum(qr.resid(q, one)^2)) <= 1e-8 * sqrt(nrow(X)))) return(X)
    g <- qr.coef(q, one)
    stay <- !is.na(g) & abs(g) > 1e-8
  } else {
    stay[ic] <- TRUE
  }
  m <- colMeans(X)
  m[stay] <- 0
  Xc <- X - rep(m, each = nrow(X))
  raw <- sqrt(colSums(X^2))
  cen <- sqrt(colSums(Xc^2))
  flat <- !stay & !(cen > design_alias_tol * raw)
  Xc[, flat] <- 0
  Xc
}",
"item 2, design_center()")
append_to("tests/testthat/test-design-offset.R",
"test_that(\"a product of shifted columns beside `0 + f` is judged centered\", {
  # `0 + f` spans the constant with no intercept column, and the rank
  # test read such a design uncentered: it dropped w1:w2 at 1e5 (capped)
  # and at 1e6 with the message \"rank deficient\", 84 to 101
  # log-likelihood units short (release review of 0.80.0)
  set.seed(7)
  n <- 300
  d <- data.frame(z1 = stats::rnorm(n), z2 = stats::rnorm(n),
                  f = factor(sample(letters[1:3], n, TRUE)))
  d$y <- c(a = 0, b = 1, c = -1)[as.character(d$f)] + 0.5 * d$z1 -
    0.4 * d$z2 + 0.6 * d$z1 * d$z2 + stats::rnorm(n, 0, 0.5)
  fz <- frm(y ~ 0 + f + z1 * z2, data = d)
  for (m in c(1e5, 1e6)) {
    d$w1 <- d$z1 + m
    d$w2 <- d$z2 + m
    cls <- character(0)
    fw <- NULL
    expect_no_message(fw <- withCallingHandlers(
      frm(y ~ 0 + f + w1 * w2, data = d),
      warning = function(x) {
        cls <<- c(cls, if (inherits(x, \"frmtmb_design_capped\")) \"capped\"
                  else class(x)[1L])
        invokeRestart(\"muffleWarning\")
      }))
    if (m == 1e5) {
      # the centered design is well conditioned here: kept and fitted
      expect_length(fw$estimates[[\"beta\"]], 6L)
      expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
                1e-6)
    } else if (length(fw$estimates[[\"beta\"]]) == 6L) {
      expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
                1e-6)
    } else {
      expect_true(\"capped\" %in% cls)
    }
  }
})",
"item 2, the reviewer's test")

# item 11: a test of merge step 3's wording of the multistart message
append_to("tests/testthat/test-nl-start.R",
"test_that(\"the multistart message names the penalized value under a prior\", {
  # merge step 3 of the 0.80.0 consolidation: what the fits maximized is
  # the penalized log-likelihood under a prior with a density, and the
  # message must say so and print it (release review of 0.80.0, item 11)
  fm <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
  pr <- set_prior(\"normal(0, 5)\", nlpar = \"b1\") +
    set_prior(\"normal(0, 5)\", nlpar = \"b2\")
  set.seed(1)
  d <- data.frame(x = runif(200, 0, 4))
  d$y <- 3 * (1 - exp(-1.2 * d$x)) + rnorm(200, 0, 0.3)
  said <- function(got) {
    m <- Filter(function(cnd) {
      inherits(cnd, c(\"frmtmb_nl_multistart\", \"frmtmb_nl_restart\"))
    }, got$conds)
    vapply(m, conditionMessage, \"\")
  }
  on <- caught(frm(fm, data = d, prior = pr))
  off <- caught(frm(fm, data = d))
  m_on <- said(on)
  m_off <- said(off)
  # whether a restart or the probe speaks depends on the optimizer's path
  skip_if(!length(m_on) || !length(m_off), \"no restart or probe message\")
  expect_true(all(grepl(\"penalized log-likelihood\", m_on, fixed = TRUE)))
  expect_false(any(grepl(\"penalized\", m_off, fixed = TRUE)))
  # the last value printed is the one the returned fit maximized
  nums <- regmatches(m_on, gregexpr(\"log-likelihood -?[0-9.]+\", m_on))
  last <- as.numeric(sub(\".*log-likelihood \", \"\", utils::tail(unlist(nums), 1)))
  lp <- as.numeric(logLik_penalized(on$value))
  expect_lt(abs(last / lp - 1), 1e-5)
})",
"item 11, the multistart wording under a prior")

# item 8: the survey's corrected-weights line
edit("dev/equiv-summary.R",
"    W2 <- vapply(strata, function(s) sum(pw[cf$status == s]), 0)",
"    # a frame newer than weights.tsv holds features it does not weigh;
    # they add no weight
    W2 <- vapply(strata, function(s) sum(pw[cf$status == s], na.rm = TRUE), 0)",
"item 8, na.rm in equiv-summary.R")

# item 7: the checklist's commands
edit("dev/release/checklist.md",
"1. Compare the library with CRAN: `Rscript dev/release/cran-lag.R` (exit 0,
   or each lagging package updated or explained).",
"1. Compare the library with CRAN:
   `Rscript dev/release/cran-lag.R \"<release library>;<user library>\"` (exit
   0, or each lagging package updated or explained).",
"item 7, step 1")
edit("dev/release/checklist.md",
"4. Run the test suites of core and of every extension, one file per process,
   with `dev/ciharden-scan.sh`:",
"4. Run the test suites of core and of every extension, one file per process,
   with `ROUND_BASE=<previous release library> bash dev/ciharden-scan.sh <run>
   ref <release library> 3`:",
"item 7, step 4")
edit("dev/release/checklist.md",
"no `LEAK` line can be untagged.",
"a `LEAK` line must have a tag.",
"item 7, step 6 wording")
edit("dev/release/checklist.md",
"10. Make sure that the installed code equals the tree: count the functions
    that differ (0 expected).",
"10. Make sure that the installed code equals the tree, with
    `dev/rel080-installed.R` (copy it and change the library): count the
    functions that differ (0 expected). Compare the top-level constants too,
    not only the functions.",
"item 7, step 10")
edit("dev/release/checklist.md",
"13. Regenerate the ported brms-suite ledger once. Report each row that flips.",
"13. Regenerate the ported brms-suite ledger once, with `dev/rel080-ledger.sh`.
    Report each row that flips.",
"item 7, step 13")

# item 11: the release commit names the fix scripts
edit("dev/rel080-commit.sh",
"resolutions), dev/rel080-release-fix.R (the brief's merge steps),
dev/rel080-news.R, dev/rel080-backlog.R.",
"resolutions), dev/rel080-release-fix.R (the brief's merge steps),
dev/rel080-release-fix2.R (conflicts found on the build, test edits),
dev/rel080-release-fix3.R (the release review's fixes),
dev/rel080-news.R, dev/rel080-backlog.R.",
"item 11, commit message")

# item 6: the backlog. Item 1 is applied, so the bernoulli mm 55 item
# stays as it is; the three missing items are added
edit("dev/test-backlog.md",
"### From the reviews, filed at the release\n\n",
paste0("### From the reviews, filed at the release\n\n",
"- **MINOR. A nonlinear fit that ends on a ridge to a linear limit can read\n",
"  finite standard errors in the body's Jacobian basis, and then nothing\n",
"  says it is not a maximum.** binomial `Vm * xp / (K + xp)` (review of\n",
"  lane nlstart, fresh row mm 22) slid from K 9.8e7 to 1.1e8 and lost the\n",
"  standard-error warning (release review of 0.80.0, item 1). 0.80.0 keeps\n",
"  the second descent off the refit of a sign probe's winner; a ridge test\n",
"  that does not go through the standard errors would close the rest.\n",
"- **MINOR. Products and powers of columns shifted by 1e6 under an ordinal\n",
"  family stop 8.6 to 28.6 short, with warnings** (release review of\n",
"  0.80.0, item 2): the thresholds take the intercept's place, so the\n",
"  centered rank test and the capped warning do not reach the design.\n",
"- **NIT. `nl_local_message()` reads `stats::logLik()`**, which stops under\n",
"  REML with a prior on the fixed effects. It is unreachable under any\n",
"  prior (release review of 0.80.0, item 9).\n"),
"item 6, the three missing backlog items")
cat("RELEASE FIX 3 DONE\n")
