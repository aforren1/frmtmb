# Fix round 2: tests for (1) newdata over the rows of the term's
# matrix, (2) both columns of the map check, (3) the rounded display of
# the default prior's ends. Part A: Rscript <this> <root>
root <- commandArgs(TRUE)[1]
edit <- function(p, pairs, append = NULL) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  if (!is.null(append)) txt <- paste0(txt, "\n", append)
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}

# (1) sar
edit("extensions/frmtmb.spatial/tests/testthat/test-sar.R", list(
  c("  expect_error(predict(fit, newdata = n3), \"with the fitted rows\",
               class = \"frmtmb_error\")
  expect_equal(unname(frm_linpred(fit, newdata = n3)),
               unname(frm_linpred(fit)[1:3]), tolerance = 1e-12)
})",
    "  expect_error(predict(fit, newdata = n3), \"one row per row of the term\",
               class = \"frmtmb_error\")
  expect_equal(unname(frm_linpred(fit, newdata = n3)),
               unname(frm_linpred(fit)[1:3]), tolerance = 1e-12)
})

test_that(\"after rows are dropped, newdata covers every row of W, as in
           brms\", {
  # brms reads newdata against W from data2, so it takes nrow(W) rows,
  # the two the fit dropped among them, and refuses the kept rows alone
  # (dev/sarfcor-fr2-rows.R)
  s <- sar_data(k = 5)
  d <- s$d
  W <- s$W
  d$x[c(4, 17)] <- NA
  fit <- suppressMessages(frm(ylag ~ x + sar(W), data = d,
                              data2 = list(W = W)))
  keep <- setdiff(seq_len(s$n), c(4, 17))
  nd <- d
  nd$x <- ifelse(is.na(nd$x), 0, nd$x) + 1
  rho <- sar_rho_of(fit)
  b <- unname(coef(fit))
  # the mean over all 25 rows of W, at the fitted rho
  want <- as.vector(solve(diag(s$n) - rho * W, b[1L] + b[2L] * nd$x))
  expect_equal(unname(fitted(fit, newdata = nd)[, \"Estimate\"]), want,
               tolerance = 1e-10)
  expect_identical(dim(predict(fit, newdata = nd)), c(as.integer(s$n), 4L))
  expect_identical(dim(as.matrix(simulate(fit, newdata = nd, nsim = 2,
                                          seed = 1))),
                   c(as.integer(s$n), 2L))
  # a missing predictor gives NA, which the lag spreads to every row,
  # as in brms
  expect_true(all(is.na(fitted(fit, newdata = transform(nd, x = d$x + 1))[
    , \"Estimate\"])))
  # the kept rows alone are refused, as brms refuses them
  expect_error(fitted(fit, newdata = nd[keep, ]),
               \"one row per row of the term's matrix: 25 rows\",
               fixed = TRUE, class = \"frmtmb_error\")
  expect_error(simulate(fit, newdata = nd[keep, ]), \"dropped 2 of those\",
               class = \"frmtmb_error\")
  # the absent case: without newdata, the fitted rows, as in brms
  expect_length(fitted(fit)[, \"Estimate\"], s$n)
})"),
  c("  # and brms refuses rows that differ, as frmtmb does
  expect_error(brms::posterior_epred(bfit, newdata = nd[1:5, ]))
})",
    "  # and brms refuses rows that differ, as frmtmb does
  expect_error(brms::posterior_epred(bfit, newdata = nd[1:5, ]))
  # after the fit drops rows: brms takes newdata over every row of W and
  # refuses the kept rows alone, and so does frmtmb
  dd <- d
  dd$x[c(5, 20)] <- NA
  keep <- setdiff(seq_len(s$n), c(5, 20))
  fd <- suppressMessages(frm(ylag ~ x + sar(Ws), data = dd,
                             data2 = list(Ws = Ws)))
  b <- unname(coef(fd))
  pd <- list(b = array(b[2L], 1),
             Intercept = b[1L] + b[2L] * mean(dd$x[keep]),
             sigma = sigma(fd), lagsar = sar_rho_of(fd))
  prior <- brms_flat_prior(bform, data = dd, family = fam,
                           data2 = list(Ws = Ws))
  code <- brms::make_stancode(bform, data = dd, family = fam,
                              prior = prior, data2 = list(Ws = Ws))
  sdat <- brms_standata(bform, data = dd, family = fam, prior = prior,
                        data2 = list(Ws = Ws))
  bd <- suppressMessages(brms::brm(bform, data = dd, family = fam,
                                   prior = prior, data2 = list(Ws = Ws),
                                   empty = TRUE))
  bd$fit <- suppressMessages(rstan::sampling(
    brms_stan_model(code), data = sdat, chains = 1, iter = 2, warmup = 0,
    algorithm = \"Fixed_param\", init = list(pd), refresh = 0, seed = 1))
  bd <- brms::rename_pars(bd)
  nf <- dd
  nf$x <- ifelse(is.na(nf$x), 0, nf$x) + 1
  be <- brms::posterior_epred(bd, newdata = nf)[1L, ]
  fe <- unname(fitted(fd, newdata = nf)[, \"Estimate\"])
  expect_length(be, s$n)
  expect_lt(max(abs(be - fe)) / max(abs(fe)), 1e-10)
  expect_error(brms::posterior_epred(bd, newdata = nf[keep, ]),
               \"must be equal to the number of observations\")
  expect_error(fitted(fd, newdata = nf[keep, ]), \"one row per row\",
               class = \"frmtmb_error\")
})")))

# (1) fcor
edit("tests/testthat/test-fcor.R", list(
  c("  expect_error(predict(fit, newdata = nd), \"with the fitted rows\",
               class = \"frmtmb_error\")",
    "  expect_error(predict(fit, newdata = nd), \"one row per row of the term\",
               class = \"frmtmb_error\")
  # after rows are dropped, newdata covers every row of M, as brms reads
  # it, and the kept rows alone are refused, as brms refuses them
  dn <- d
  dn$x[c(3, 7)] <- NA
  fn <- suppressMessages(frm(y ~ x + fcor(V), data = dn,
                             data2 = list(V = V)))
  nf <- dn
  nf$x <- ifelse(is.na(nf$x), 0, nf$x) + 1
  expect_equal(unname(fitted(fn, newdata = nf)[, \"Estimate\"]),
               unname(frm_linpred(fn, newdata = nf)), tolerance = 1e-12)
  expect_length(fitted(fn, newdata = nf)[, \"Estimate\"], fd$n)
  sn2 <- as.matrix(simulate(fn, newdata = nf, nsim = 2, seed = 1))
  expect_identical(dim(sn2), c(as.integer(fd$n), 2L))
  expect_error(fitted(fn, newdata = nf[-c(3, 7), ]),
               \"one row per row of the term's matrix: 10 rows\",
               fixed = TRUE, class = \"frmtmb_error\")")))

# (2) the map check reads both columns
edit("tests/testthat/test-autocor-seam.R", list(
  c("  # whiten() whose log-Jacobian is one value per row",
    "  # unwhiten() that is right on one column of innovations and wrong on
  # the others
  reg(\"col1un\", unwhiten = function(u, mu, sigma, theta, aux) {
    out <- mu + sigma * (aux$L %*% u)
    if (is.matrix(u) && ncol(u) > 1L) out[, -1L] <- 0
    out
  })
  # whiten() whose log-Jacobian is one value per row"),
  c("  expect_error(fit_with(\"vecjac\"), \"must give one finite log-Jacobian\",
               fixed = TRUE, class = \"frmtmb_error\")",
    "  expect_error(fit_with(\"vecjac\"), \"must give one finite log-Jacobian\",
               fixed = TRUE, class = \"frmtmb_error\")
  expect_error(fit_with(\"col1un\"),
               \"unwhiten() of a matrix of innovations must give\",
               fixed = TRUE, class = \"frmtmb_error\")")))

# (3) the announced ends are rounded, the prior's are not
edit("extensions/frmtmb.sample/tests/testthat/test-autocor-seam-draws.R",
     list(), append = "
test_that(\"the announced interval of a default prior is rounded, and the
           prior itself is not\", {
  frmtmb::frmtmb_local_registries()
  reg_lagterm()
  cs <- lag_case()
  fit <- frm(y ~ x + lagtest(W), data = cs$d, data2 = list(W = cs$W))
  ap <- frmtmb::frmtmb_autocor_pars(fit)
  dp <- frmtmb.sample:::default_priors_for(fit)
  s <- Filter(function(e) identical(e$dpar, \"lagsar\"), unclass(dp))[[1L]]
  # the prior holds the interval exactly
  expect_identical(c(s$dist$lower, s$dist$upper), c(ap$lower, ap$upper))
  msgs <- character()
  withCallingHandlers(
    frmtmb.sample:::announce_default_priors(dp, NULL),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart(\"muffleMessage\")
    })
  ln <- grep(\"lagsar \", unlist(strsplit(msgs, \"\\n\", fixed = TRUE)),
             value = TRUE, fixed = TRUE)
  expect_length(ln, 1L)
  shown <- regmatches(ln, regexec(\"uniform[(]([^,]+), ([^)]+)[)]\", ln))[[1L]]
  expect_length(shown, 3L)
  # six digits at most, and each within a relative 1e-5 of the end
  expect_true(all(nchar(gsub(\"[^0-9]\", \"\", shown[2:3])) <= 6L))
  ends <- as.numeric(shown[2:3])
  expect_lt(max(abs(ends - c(ap$lower, ap$upper)) /
                  abs(c(ap$lower, ap$upper))), 1e-5)
})")
cat("edited", root, "\n")
