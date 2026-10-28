# Replace one test_that() block of a test file, found by its title, with
# the text of another file. Usage inside R: replace_test(file, title, new)
replace_test <- function(f, title, new) {
  x <- readLines(f)
  s <- which(startsWith(x, paste0("test_that(\"", title)))
  stopifnot(length(s) == 1L)
  e <- s + which(x[(s + 1L):length(x)] == "})")[1L]
  writeLines(c(x[seq_len(s - 1L)], new, x[(e + 1L):length(x)]), f)
}

f <- "extensions/frmtmb.sample/tests/testthat/test-default-priors-autocor.R"
replace_test(f, "the default on sigma is brms's, on the innovation sd", c(
"test_that(\"the default on sigma is brms's, on sigma itself\", {",
"  # brms's cholesky_cor_ar1(), _ma1() and _arma1() are covariances, so",
"  # its sigma under a cov = TRUE term is the innovation sd, and so is",
"  # frmtmb's since sigma became the innovation sd. The default",
"  # student_t(3, 0, s) is brms's density at sigma = exp(eta), with the",
"  # log link's Jacobian and no factor of the correlation parameters.",
"  # frmtmb.sample 0.23.0 read it at sigma * c with the Jacobian of c,",
"  # because frmtmb's sigma was the marginal sd.",
"  d <- dpa_data()",
"  for (fo in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + ma(t, g, cov = TRUE),",
"                  y ~ x + arma(t, g, cov = TRUE),",
"                  y ~ x + ar(t, g, p = 2, cov = TRUE), y ~ x + cosy(t, g))) {",
"    fit <- frm(fo, data = d, dry_run = \"objective\")",
"    ri <- dpa_ri(fit)",
"    ac <- fit$frame$autocor[[1L]]",
"    sg <- Filter(function(e) !identical(e$comp, \"thetaac\") &&",
"                   identical(e$dist$kind, \"t\") && identical(e$scale, \"sd\"),",
"                 ri$entries)",
"    expect_length(sg, 1L)",
"    expect_null(sg[[1L]]$innov)",
"    sc <- sg[[1L]]$dist$scale",
"    nlp <- neg_log_prior_fn(sg)",
"    set.seed(2)",
"    v <- vapply(1:5, function(k) {",
"      p <- fit$frame$par_template",
"      p[[\"thetaac\"]][] <- stats::rnorm(ac$npar)",
"      eta <- stats::rnorm(1, -0.5, 0.5)",
"      p[[sg[[1L]]$comp]][sg[[1L]]$idx] <- eta",
"      # brms's density on sigma = exp(eta), carried to eta",
"      ref <- stats::dt(exp(eta) / sc, 3, log = TRUE) - log(sc) + eta",
"      c(nlp(p), ref)",
"    }, c(0, 0))",
"    expect_lt(diff(range(colSums(v))) / max(abs(v)), 1e-10,",
"              label = deparse(fo))",
"  }",
"  # under sigma ~ x the default is on brms's Intercept_sigma, the log",
"  # innovation sd at the predictor means, which no correlation",
"  # parameter moves",
"  fx <- frm(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), data = d,",
"            dry_run = \"objective\")",
"  ex <- Filter(function(e) identical(e$comp, \"betad\"), dpa_ri(fx)$entries)",
"  expect_length(ex, 1L)",
"  nlx <- neg_log_prior_fn(ex)",
"  set.seed(4)",
"  p <- fx$frame$par_template",
"  p[[\"betad\"]][] <- c(-0.4, 0.3)",
"  vx <- vapply(1:5, function(k) {",
"    p[[\"thetaac\"]][] <- stats::rnorm(1, 0, 2)",
"    nlx(p)",
"  }, 0)",
"  expect_identical(diff(range(vx)), 0)",
"  msg <- dpa_messages(frmtmb.sample:::sample_resolve_priors(",
"    frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = \"objective\"),",
"    NULL))",
"  expect_false(any(grepl(\"innovation sd\", msg, fixed = TRUE)))",
"})"))
