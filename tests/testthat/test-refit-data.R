# Refits read the data the fit holds, not the model frame and not the
# workspace.
#
# drop1(), update() without data, influence() and step() evaluate the
# formula again. The model frame holds `log(x2)`, not `x2`, so up to
# 0.71.0 drop1() stopped with "The model uses `x2`, which is not a
# column of `data`" on any refit that kept a transformed term, and
# update() read whatever object the call's `data` symbol named at the
# time. dev/drop1raw-findings.md has the survey, the comparison with
# stats::drop1() and lme4 and the measurements of what the fit keeps.

refit_data_case <- function(n = 200, seed = 5) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                  g = factor(sample(letters[1:8], n, TRUE)))
  d$z <- 1 + d$x + log(d$x2) + 0.3 * d$x * log(d$x2) + rnorm(8)[d$g] +
    rnorm(n)
  d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x + 0.3 * log(d$x2)))
  d
}

# frmtmb's LRT and Df per dropped term, from a fit on `data`. The call
# holds its arguments by value, so the only thing a refit has to find is
# what the formula reads
refit_drop1 <- function(f, data, family = gaussian(), scope = NULL) {
  fit <- suppressMessages(do.call(frm, list(formula = bf(f),
                                            family = family, data = data)))
  t <- if (is.null(scope)) {
    drop1(fit, test = "Chisq")
  } else {
    drop1(fit, scope = scope, test = "Chisq")
  }
  list(rows = rownames(t)[-1L], lrt = t$LRT[-1L], df = t$Df[-1L])
}

# stats::drop1() on lm() drops the fit's own design columns; with
# scale = 0 its likelihood-ratio statistic is n log(RSS_i / RSS_0)
lm_drop1 <- function(f, data, scope = NULL) {
  m <- stats::lm(f, data = data)
  t <- if (is.null(scope)) stats::drop1(m) else
    stats::drop1(m, scope = scope)
  list(rows = rownames(t)[-1L], df = t$Df[-1L],
       lrt = stats::nobs(m) * log(t$RSS[-1L] / t$RSS[1L]))
}

# the statistics agree to a part in 1e6 of their own size; measured at
# 1e-12 or better (dev/drop1raw-compare.R)
expect_same_drop1 <- function(a, b) {
  expect_identical(a$rows, b$rows)
  expect_identical(as.numeric(a$df), as.numeric(b$df))
  expect_lt(max(abs(a$lrt - b$lrt) / abs(b$lrt)), 1e-6)
}

test_that("drop1() refits terms that read a raw variable, as lm() does", {
  d <- refit_data_case()
  fs <- list(z ~ x + log(x2), z ~ scale(x) + x2, z ~ I(x^2) + x2,
             z ~ x + splines::ns(x2, 3), z ~ poly(x2, 2) + x,
             z ~ x * log(x2), log(z + 10) ~ x + log(x2))
  for (f in fs) {
    expect_same_drop1(refit_drop1(f, d), lm_drop1(f, d))
  }
})

test_that("drop1(scope = ) refits the terms it keeps", {
  d <- refit_data_case()
  f <- z ~ x + log(x2) + x3
  a <- refit_drop1(f, d, scope = ~ x)
  expect_identical(a$rows, "x")
  expect_same_drop1(a, lm_drop1(f, d, scope = ~ x))
})

test_that("drop1() on a poisson model matches glm()'s deviance", {
  d <- refit_data_case()
  f <- cnt ~ scale(x) + splines::ns(x2, 3)
  g <- stats::glm(f, family = poisson(), data = d)
  t <- stats::drop1(g, test = "Chisq")
  expect_same_drop1(refit_drop1(f, d, poisson()),
                    list(rows = rownames(t)[-1L], df = t$Df[-1L],
                         lrt = t$LRT[-1L]))
})

test_that("drop1() on a mixed model matches lme4", {
  skip_if_not_installed("lme4")
  d <- refit_data_case()
  f <- z ~ scale(x) + splines::ns(x2, 3) + (1 | g)
  # lme4 refits through update(), which evaluates the call's `data`
  # symbol in the formula's environment
  e <- new.env(parent = globalenv())
  e$d <- d
  m <- eval(quote(lme4::lmer(z ~ scale(x) + splines::ns(x2, 3) + (1 | g),
                             data = d, REML = FALSE)), e)
  t <- stats::drop1(m, test = "Chisq")
  expect_same_drop1(refit_drop1(f, d),
                    list(rows = rownames(t)[-1L], df = t$npar[-1L],
                         lrt = t$LRT[-1L]))
})

test_that("drop1() keeps the fit's rows and bases when na.action drops rows", {
  # ns(), scale() and poly() were evaluated on every row before the
  # missing x3 removed three; each submodel keeps those bases, so it is
  # nested in the fit, as stats::drop1() keeps the fit's design columns
  d <- refit_data_case()
  d$x3[c(3, 17, 40)] <- NA
  for (f in list(z ~ splines::ns(x2, 3) + x3, z ~ scale(x) + x3,
                 z ~ poly(x2, 2) + x3)) {
    expect_same_drop1(refit_drop1(f, d), lm_drop1(f, d))
  }
})

test_that("drop1() reads an interval cens() bound the frame leaves out", {
  d <- refit_data_case()
  d$cc <- ifelse(seq_len(nrow(d)) %% 5 == 0, "interval", "none")
  d$z2 <- d$z + 1
  fit <- frm(bf(z | cens(cc, z2) ~ x + x3), data = d)
  t <- drop1(fit)
  expect_identical(rownames(t), c("<none>", "x", "x3"))
  expect_false(anyNA(t$AIC))
})

test_that("drop1() needs nothing from where the fit was made", {
  d <- refit_data_case()
  # the family and the data are a function's locals, the formula is not
  make <- function(form, dat) {
    fam <- poisson()
    frm(bf(form), family = fam, data = dat)
  }
  fit <- make(cnt ~ x + log(x2), d)
  expect_identical(rownames(drop1(fit)), c("<none>", "x", "log(x2)"))
  # a bound on the dropped coefficient goes with it, and the message
  # says so; a bound adds no density, so this is still a likelihood
  # ratio
  fp <- frm(bf(z ~ x + x3), data = d,
            prior = set_prior("", class = "b", coef = "x3", lb = -10))
  expect_message(t <- drop1(fp), "drop1[(][)] drops the prior")
  expect_identical(rownames(t), c("<none>", "x", "x3"))
})

test_that("a fit whose data are gone refits as one whose data are not", {
  d <- refit_data_case()
  seen <- frm(bf(z ~ x + log(x2) + x3), data = d)
  gone <- local({
    d_gone <- d
    frm(bf(z ~ x + log(x2) + x3), data = d_gone)
  })
  expect_false(exists("d_gone"))
  expect_identical(drop1(gone)$AIC, drop1(seen)$AIC)
  u_seen <- update(seen, . ~ . - x3)
  u_gone <- update(gone, . ~ . - x3)
  expect_identical(as.numeric(logLik(u_gone)), as.numeric(logLik(u_seen)))
  # the new fit's call names the data as the original call did
  expect_identical(u_gone$call$data, as.name("d_gone"))
})

test_that("update() refits on the fit's data after the object changes", {
  # brms's update() reads the data the fit holds; with x2 changed in the
  # workspace, the refit is still the fitted model
  d <- refit_data_case()
  fit <- frm(bf(z ~ x + log(x2)), data = d)
  d$x2 <- d$x2^2 + 1
  d <- d[-(1:10), ]
  u <- update(fit)
  expect_identical(nobs(u), nobs(fit))
  expect_identical(as.numeric(logLik(u)), as.numeric(logLik(fit)))
  # a variable the fit does not hold is not read from changed data,
  # which would read every other variable again too: brms's refusal
  d$x4 <- rnorm(nrow(d))
  expect_error(update(fit, . ~ . + x4),
               paste0("New variables found: 'x4'\nPlease supply your data ",
                      "again via argument 'newdata'."), fixed = TRUE)
  # and `newdata` still replaces the data
  expect_identical(nobs(update(fit, newdata = d[1:50, ])), 50L)
  expect_identical(nobs(update(fit, . ~ . + x4, newdata = d)), nrow(d))
})

test_that("a new variable is read from data that still hold the fit's", {
  # stats::update() and add1() read it from the call's data; that is
  # done when those data are the fit's, rows and values
  d <- refit_data_case()
  d$x3[c(3, 17)] <- NA
  fit <- suppressMessages(frm(bf(z ~ x + log(x2) + x3), data = d))
  d$x4 <- rnorm(nrow(d))
  u4 <- suppressMessages(update(fit, . ~ . + x4))
  ref <- suppressMessages(frm(bf(z ~ x + log(x2) + x3 + x4), data = d))
  expect_identical(nobs(u4), nobs(fit))
  expect_identical(as.numeric(logLik(u4)), as.numeric(logLik(ref)))
  a1 <- suppressMessages(stats::add1(fit, ~ . + x4 + I(x^2)))
  expect_identical(rownames(a1), c("<none>", "x4", "I(x^2)"))
  # one value changed is another model
  d$x[5] <- d$x[5] + 1
  expect_error(update(fit, . ~ . + x4), "New variables found: 'x4'",
               fixed = TRUE)
  expect_error(stats::add1(fit, ~ . + x4), "New variables found",
               fixed = TRUE)
})

test_that("update() keeps the rows a univariate subset() leaves out", {
  # brms's nobs() counts them, and so did update() when it read the
  # call's data; the fit's own data hold them for the same reason
  set.seed(11)
  d <- data.frame(x = rnorm(160), x2 = rexp(160) + 0.1,
                  sb = rep(c(TRUE, TRUE, FALSE), length.out = 160))
  d$z <- 1 + d$x + log(d$x2) + rnorm(160)
  # a missing value on a row the subset leaves out costs no row
  d$x[c(3, 6)] <- NA
  fit <- suppressMessages(frm(bf(z | subset(sb) ~ x + log(x2)), data = d))
  u <- local({
    d <- NULL
    suppressMessages(update(fit))
  })
  expect_identical(nobs(u), nobs(fit))
  expect_identical(nobs(fit), 160L)
  expect_identical(u$frame$n_obs, fit$frame$n_obs)
  expect_identical(as.numeric(logLik(u)), as.numeric(logLik(fit)))
})

test_that("a saved fit refits in a new session that has no data", {
  skip_on_cran()
  d <- refit_data_case()
  fit <- frm(bf(z ~ x + log(x2) + x3), data = d)
  rds <- tempfile(fileext = ".rds")
  res <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(rds, res, script)), add = TRUE)
  saveRDS(fit, rds)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    "suppressPackageStartupMessages(library(frmtmb))",
    sprintf("f <- readRDS(%s)", deparse(rds)),
    "st <- stats::step(f, trace = 0)",
    "out <- list(has_d = exists('d'), aic = drop1(f)$AIC,",
    "  ll = as.numeric(logLik(stats::update(f, . ~ . - x3))),",
    "  step = deparse1(formula(st)), step_ll = as.numeric(logLik(st)))",
    sprintf("saveRDS(out, %s)", deparse(res))), script)
  log <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                  c("--vanilla", shQuote(script)),
                                  stdout = TRUE, stderr = TRUE))
  expect_true(file.exists(res), info = paste(log, collapse = "\n"))
  out <- readRDS(res)
  expect_false(out$has_d)
  # a fresh process may round differently in the last bits (rtmb-pitfalls
  # item 21), so the comparison is to a part in 1e8
  ll <- as.numeric(logLik(update(fit, . ~ . - x3)))
  expect_lt(max(abs(out$aic - drop1(fit)$AIC) / abs(drop1(fit)$AIC)), 1e-8)
  expect_lt(abs(out$ll - ll) / abs(ll), 1e-8)
  st <- stats::step(fit, trace = 0)
  expect_identical(out$step, deparse1(formula(st)))
  expect_lt(abs(out$step_ll - as.numeric(logLik(st))) /
              abs(as.numeric(logLik(st))), 1e-8)
})

test_that("step() removes a term", {
  # step() passes its update as text, "~ . - x3", which update() used to
  # leave unevaluated, so step() stopped with "object 'change' not found"
  set.seed(5)
  d <- data.frame(x = rnorm(120), x3 = rnorm(120))
  d$z <- 1 + d$x + rnorm(120)
  fit <- frm(bf(z ~ x + x3), data = d)
  st <- stats::step(fit, trace = 0)
  expect_identical(deparse1(formula(st)), "z ~ x")
  expect_identical(deparse1(formula(stats::step(stats::lm(z ~ x + x3, d),
                                                trace = 0))),
                   "z ~ x")
})

step_case <- function(seed = 21, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x4 = rnorm(n))
  d$z <- 1 + d$x + log(d$x2) + rnorm(n)
  d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x + 0.3 * log(d$x2)))
  d
}

test_that("step() moves on the fit's data, not on the workspace's", {
  # stats::step() evaluates the call update(evaluate = FALSE) returns;
  # that call named `d`, so a changed `d` became the model and a removed
  # one stopped step() after drop1() had chosen the move
  d <- step_case()
  fit <- frm(bf(z ~ x + log(x2) + x4), data = d)
  ref <- stats::step(fit, trace = 0)
  # the same path as lm's
  expect_identical(deparse1(formula(ref)),
                   deparse1(formula(stats::step(
                     stats::lm(z ~ x + log(x2) + x4, d), trace = 0))))
  expect_false(identical(deparse1(formula(ref)), "z ~ x + log(x2) + x4"))
  ll <- as.numeric(logLik(ref))
  d$x2 <- d$x2^2 + 1
  st <- stats::step(fit, trace = 0)
  expect_identical(deparse1(formula(st)), deparse1(formula(ref)))
  expect_identical(as.numeric(logLik(st)), ll)
  rm(d)
  st <- stats::step(fit, trace = 0)
  expect_identical(as.numeric(logLik(st)), ll)
  # the chosen model's call names the data as the fit's call did
  expect_identical(st$call$data, as.name("d"))
})

test_that("step() on a poisson fit takes glm()'s path", {
  d <- step_case()
  fit <- frm(bf(cnt ~ x + log(x2) + x4), family = poisson(), data = d)
  g <- stats::step(stats::glm(cnt ~ x + log(x2) + x4, family = poisson(),
                              data = d), trace = 0)
  st <- stats::step(fit, trace = 0)
  expect_identical(deparse1(formula(st)), deparse1(formula(g)))
  # glm's deviance and frmtmb's log-likelihood differ by a constant of
  # the data, which cancels between the two fits
  full <- frm(bf(cnt ~ x + log(x2) + x4), family = poisson(), data = d)
  dll <- 2 * (as.numeric(logLik(full)) - as.numeric(logLik(st)))
  ddev <- stats::deviance(g) -
    stats::deviance(stats::glm(cnt ~ x + log(x2) + x4,
                               family = poisson(), data = d))
  expect_lt(abs(dll - ddev) / max(abs(ddev), 1e-3), 1e-6)
})

test_that("update() and drop1() keep the fit's bases when rows are dropped", {
  # the fit evaluates ns() on every row before na.action removes three;
  # a refit on the fit's rows has to keep those knots, or update(fit)
  # is another model and drop1() and update() disagree on a submodel
  d <- refit_data_case()
  d$x3[c(3, 17, 40)] <- NA
  fit <- suppressMessages(frm(bf(z ~ splines::ns(x2, 3) + x3), data = d))
  ll <- as.numeric(logLik(fit))
  # measured: 1.49e-2 apart when the refit computed the knots again
  expect_lt(abs(as.numeric(logLik(update(fit))) - ll) / abs(ll), 1e-10)
  u <- update(fit, ~ . - x3)
  a <- drop1(fit)["x3", "AIC"]
  expect_lt(abs(a - stats::AIC(u)) / stats::AIC(u), 1e-10)
  # poly() and scale(): the same coefficients, under the same names
  for (f in list(z ~ poly(x2, 2) + x3, z ~ scale(x) + x3)) {
    fp <- suppressMessages(frm(bf(f), data = d))
    b0 <- fixef(fp)[, "Estimate"]
    b1 <- fixef(update(fp))[, "Estimate"]
    expect_identical(names(b1), names(b0))
    expect_lt(max(abs(b1 - b0) / pmax(abs(b0), 1)), 1e-6)
  }
})

test_that("drop1() keeps a kept term's name, so its prior stays", {
  # a bound on b_scalex: the refit without x3 keeps scale(x), and its
  # coefficient is still scalex, so that refit keeps the bound
  d <- refit_data_case()
  fs <- frm(bf(z ~ scale(x) + x3), data = d,
            prior = set_prior("", class = "b", coef = "scalex", lb = -100))
  msgs <- character(0)
  t <- withCallingHandlers(drop1(fs), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  # one message: the refit without scale(x) has no b_scalex to bound
  expect_identical(sum(grepl("b_scalex", msgs, fixed = TRUE)), 1L)
  expect_identical(rownames(t), c("<none>", "scale(x)", "x3"))
})

test_that("update() takes the formula, family and control from the fit", {
  # brms takes them from the object; the call's symbols name objects of
  # the session that made the fit
  d <- refit_data_case()
  fit_with <- function(form, fam, data) frm(form, family = fam, data = data)
  fit <- fit_with(bf(z ~ x + x3), gaussian(), d)
  ll <- as.numeric(logLik(fit))
  expect_identical(as.numeric(logLik(update(fit))), ll)
  u <- update(fit, ~ . - x3)
  expect_identical(deparse1(formula(u)), "z ~ x")
  # objects of the same names here do not become the model
  form <- bf(z ~ x)
  fam <- student()
  u2 <- update(fit)
  expect_identical(deparse1(formula(u2)), "z ~ x + x3")
  expect_identical(family(u2)[["family"]], "gaussian")
  expect_identical(as.numeric(logLik(u2)), ll)
  # a family given to update() replaces the fit's
  expect_identical(family(update(fit, family = student()))[["family"]],
                   "student")
  # and the call reads as the original one
  expect_identical(u2$call$formula, as.name("form"))
})

test_that("drop1() and anova() refuse a fit whose prior adds a density", {
  # logLik() of such a fit is log-likelihood plus log prior, so the
  # difference is not a likelihood ratio and can be negative
  d <- refit_data_case()
  pr <- set_prior("normal(0, 1)", class = "b")
  fp <- frm(bf(z ~ x + x3), data = d, prior = pr)
  expect_error(drop1(fp), "is not a likelihood ratio", fixed = TRUE)
  expect_error(stats::step(fp, trace = 0), "is not a likelihood ratio",
               fixed = TRUE)
  f0 <- frm(bf(z ~ x), data = d)
  expect_error(anova(f0, fp), "model 2 has a prior", fixed = TRUE)
  # a prior that only sets bounds adds nothing: logLik() is the ML's
  fb <- frm(bf(z ~ x + x3), data = d,
            prior = set_prior("", class = "b", coef = "x3", lb = -10))
  fm <- frm(bf(z ~ x + x3), data = d)
  expect_lt(abs(as.numeric(logLik(fb)) - as.numeric(logLik(fm))) /
              abs(as.numeric(logLik(fm))), 1e-8)
  expect_message(t <- drop1(fb), "drops the prior specification")
  expect_identical(rownames(t), c("<none>", "x", "x3"))
  expect_identical(nrow(anova(f0, fb)), 2L)
})

test_that("update() on an na.exclude fit pads as the fit does", {
  # lm's update() keeps 200 rows of residuals; the refit on the fit's
  # rows takes the fit's missing rows
  d <- refit_data_case()
  d$x3[c(3, 17, 40)] <- NA
  fit <- suppressMessages(frm(bf(z ~ x + x3), data = d,
                              na.action = stats::na.exclude))
  u <- update(fit, ~ . - x3)
  l <- stats::update(stats::lm(z ~ x + x3, data = d,
                               na.action = stats::na.exclude), ~ . - x3)
  expect_identical(nrow(residuals(fit)), 200L)
  expect_identical(nrow(residuals(u)), length(residuals(l)))
  expect_identical(which(is.na(residuals(u)[, 1])), c(3L, 17L, 40L))
})

test_that("an edited update() call records the model it fits", {
  # the call update(evaluate = FALSE) returns names the refit it was
  # made for; edited, it is another model, and the old fit's call,
  # bases and missing rows do not belong to it
  d <- refit_data_case()
  d$x3[c(3, 17, 40)] <- NA
  fit <- suppressMessages(frm(bf(z ~ splines::ns(x2, 3) + x + x3),
                              data = d))
  cl <- update(fit, evaluate = FALSE)
  cl$formula <- bf(z ~ x)
  e <- suppressMessages(eval(cl))
  expect_identical(deparse1(formula(e)), "z ~ x")
  expect_false(identical(getCall(e)$formula, getCall(fit)$formula))
  # the refit data hold 197 complete rows, and none was dropped
  expect_null(stats::na.action(e))
  # unedited, the call is the update
  u <- suppressMessages(eval(update(fit, evaluate = FALSE)))
  expect_identical(as.numeric(logLik(u)), as.numeric(logLik(fit)))
  expect_identical(getCall(u)$formula, getCall(fit)$formula)
})

test_that("data taken from an update() call fit as plain data", {
  d <- refit_data_case()
  d$x3[c(3, 17, 40)] <- NA
  fit <- suppressMessages(frm(bf(z ~ splines::ns(x2, 3) + x + x3),
                              data = d))
  dd <- update(fit, evaluate = FALSE)$data
  plain <- dd
  attr(plain, "frmtmb_refit") <- NULL
  # the ns() knots come from these 197 rows, not from the old fit's 200
  f1 <- frm(bf(z ~ splines::ns(x2, 3)), data = dd)
  f0 <- frm(bf(z ~ splines::ns(x2, 3)), data = plain)
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  expect_null(stats::na.action(f1))
  expect_identical(getCall(f1)$data, quote(dd))
})

test_that("data that lack a variable the fit read are not its data", {
  # a column missing from the call's data would be read from wherever
  # R finds it: here a vector of the same name and length
  d <- refit_data_case()
  fit <- frm(bf(z ~ x + log(x2)), data = d)
  d$w <- rnorm(nrow(d))
  x2 <- d$x2 * 3 + 1
  d$x2 <- NULL
  expect_error(update(fit, . ~ . + w), "New variables found: 'w'",
               fixed = TRUE)
  # with the column back the data are the fit's again
  d$x2 <- fit$frame$raw_vars$x2
  expect_identical(nobs(update(fit, . ~ . + w)), nobs(fit))
})

test_that("influence() keeps the fit's bases, as lme4 does", {
  skip_if_not_installed("lme4")
  d <- refit_data_case()
  # lme4 refits with update(subset = ), and model.frame() evaluates
  # poly() on every row before it subsets
  e <- new.env(parent = globalenv())
  e$d <- d
  m <- eval(quote(lme4::lmer(z ~ poly(x2, 2) + x + (1 | g), data = d,
                             REML = FALSE)), e)
  ref <- stats::influence(m, groups = "g")[["fixed.effects[-g]"]]
  fit <- frm(bf(z ~ poly(x2, 2) + x + (1 | g)), data = d)
  dfb <- max(abs(sweep(ref, 2, lme4::fixef(m))))
  # a part in 1e4 of the largest displacement; measured at 4.6e-6, and
  # 1.3 when each deletion computed the basis again
  # (dev/drop1raw-influence.R)
  close <- function(inf) {
    fx <- inf$fixed[rownames(ref), seq_len(ncol(ref))]
    max(abs(fx - ref)) / dfb
  }
  # the data given, as the original data
  expect_lt(close(influence(fit, groups = "g", data = d)), 1e-4)
  # the default data is the fit's own, with x2 for poly() to read
  expect_lt(close(influence(fit, groups = "g")), 1e-4)
})

test_that("a re_formula term and make_conditions() read the raw variables", {
  d <- refit_data_case()
  fit <- suppressMessages(frm(bf(z ~ x + log(x2) + (1 + log(x2) | g)),
                              data = d))
  # the fit's own term: its columns were read off the model frame, where
  # log(x2) has no x2 to be evaluated from, and the term matched nothing
  expect_identical(fitted(fit, re_formula = ~ (1 + log(x2) | g)),
                   fitted(fit))
  # brms's make_conditions() reads the stored data, which holds x2
  mc <- make_conditions(fit, "x2")
  expect_identical(mc$x2, make_conditions(d, "x2")$x2)
})

test_that("a list keeps the raw variables a refit reads", {
  d <- refit_data_case(n = 60)
  # a list has no row names to match the frame's against
  fl <- frm(bf(z ~ x + log(x2)), data = as.list(d))
  expect_identical(rownames(drop1(fl)), c("<none>", "x", "log(x2)"))
  expect_identical(drop1(fl)$AIC,
                   drop1(frm(bf(z ~ x + log(x2)), data = d))$AIC)
})

test_that("a tibble with a transformed term fits without a warning", {
  skip_if_not_installed("tibble")
  d <- refit_data_case(n = 60)
  # setting row names on a tibble warned at every fit with a transform
  expect_no_warning(ft <- frm(bf(z ~ x + log(x2)),
                              data = tibble::as_tibble(d)))
  expect_identical(drop1(ft)$AIC,
                   drop1(frm(bf(z ~ x + log(x2)), data = d))$AIC)
})
