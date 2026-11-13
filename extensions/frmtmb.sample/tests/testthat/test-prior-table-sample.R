# frm_sample()'s reading of prior tables and its default-prior message,
# lane priortab of the round of 2026-11-11. brms 2.23.0's
# default_prior() bounds sigma below by 0 and student()'s nu by 1, and
# lists a shifted lognormal's ndt on (0, min_Y); through frmtmb 0.81.0
# the table's lb and ub were NA there.

pts_data <- function() {
  set.seed(4)
  n <- 60
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(1:6, each = 10)))
  d$y <- 1 + 0.4 * d$z + stats::rnorm(6, 0, 0.4)[d$g] + stats::rnorm(n)
  d$yp <- exp(0.2 * d$z + stats::rnorm(n, 0, 0.3)) + 0.3
  d
}

pts_rp <- function(fit, prior) {
  sample_resolve_priors(fit, prior)
}

test_that("the sampling table's lb and ub are brms's", {
  d <- pts_data()
  tab <- as.data.frame(default_prior(y ~ z + (1 | g), data = d,
                                     family = student(), route = "sample"))
  expect_identical(tab$lb[tab$class == "sigma"], 0)
  expect_identical(tab$lb[tab$class == "nu"], 1)
  expect_identical(tab$lb[tab$class == "sd" & tab$group == ""], 0)
  sl <- as.data.frame(default_prior(yp ~ z, data = d,
                                    family = shifted_lognormal(),
                                    route = "sample"))
  nd <- sl[sl$class == "ndt", ]
  expect_identical(c(nd$lb, nd$ub), c(0, min(d$yp)))
})

test_that("a coerced sampling table edited to brms's (flat) is read", {
  d <- pts_data()
  f <- y ~ z + (1 | g)
  fit <- frm(f, data = d, dry_run = "objective")
  ps <- as.data.frame(default_prior(f, data = d, route = "sample"))
  msg <- "that as.data.frame() turned into a plain data frame"
  # brms's own ways to ask for (flat): "" in the cell, or set_prior("")
  ed <- ps
  ed$prior[ed$class == "sd" & ed$group == ""] <- ""
  for (pr in list(as.brmsprior(ed),
                  as.brmsprior(ps[ps$class != "sigma", ]) +
                    set_prior("", class = "sigma"))) {
    rp <- suppressMessages(pts_rp(fit, pr))
    expect_false(is.null(rp$ri))
  }
  # the guard's present case: a coerced route = "fit" table, and a frame
  # written by hand with frmtmb's printed (flat) on a slot brms gives a
  # default, are still refused
  pf <- as.data.frame(default_prior(f, data = d))
  for (pr in list(as.brmsprior(pf),
                  as.brmsprior(data.frame(prior = "(flat)",
                                          class = "sigma")))) {
    expect_error(suppressMessages(pts_rp(fit, pr)), msg, fixed = TRUE)
  }
  # the decision is per row (the review's two frames): a coerced fit
  # table with one row set to brms's own default still has its other
  # (flat) rows refused, and so does a coerced sampling table with
  # frmtmb's printed "(flat)" written into a row that had a default.
  # Each refused row is named
  h1 <- pf
  h1$prior[h1$class == "sigma"] <- "student_t(3, 0, 2.5)"
  e1 <- tryCatch(suppressMessages(pts_rp(fit, as.brmsprior(h1))),
                 error = conditionMessage)
  expect_match(e1, msg, fixed = TRUE)
  expect_match(e1, "class=Intercept", fixed = TRUE)
  expect_match(e1, "class=sd", fixed = TRUE)
  expect_false(grepl("class=sigma", e1, fixed = TRUE))
  h4 <- ps
  h4$prior[h4$class == "Intercept"] <- "(flat)"
  e4 <- tryCatch(suppressMessages(pts_rp(fit, as.brmsprior(h4))),
                 error = conditionMessage)
  expect_match(e4, "(class=Intercept)", fixed = TRUE)
  expect_false(grepl("class=sd", e4, fixed = TRUE))
})

test_that("frm_sample(fit) reads a fit's bounds-only prior as given", {
  # a bounds-only specification the user wrote is kept on the fit, so
  # the sampling run from the fit and the one given the prior directly
  # resolve one prior: brms's "" (flat) on that slot
  d <- pts_data()
  ml <- frm(y ~ z + (1 | g), data = d)
  rows <- function(rp) {
    x <- as.data.frame(rp$effective)
    x[order(x$class), c("prior", "class", "lb")]
  }
  for (cl in c("sd", "sigma")) {
    sp <- set_prior("", class = cl, lb = 0)
    f <- frm(y ~ z + (1 | g), data = d, prior = sp)
    via_fit <- suppressMessages(sample_resolve_priors(f, NULL,
                                                      base = f$prior))
    direct <- suppressMessages(pts_rp(ml, sp))
    expect_identical(rows(via_fit), rows(direct))
    expect_identical(rows(via_fit)$prior[rows(via_fit)$class == cl], "")
  }
})

test_that("a table listing Ya for Y_a replaces the defaults of its slots", {
  d <- pts_data()
  d$Y_a <- d$y
  d$Y_b <- d$y + stats::rnorm(nrow(d))
  fm <- bf(mvbind(Y_a, Y_b) ~ z) + set_rescor(FALSE)
  fit <- frm(fm, data = d, dry_run = "objective")
  tab <- default_prior(fm, data = d, route = "sample")
  eff <- as.data.frame(suppressMessages(pts_rp(fit, tab))$effective)
  # one row per slot: the table's 8, no default listed beside them
  expect_identical(nrow(eff), nrow(tab))
  expect_setequal(unique(eff$resp), c("Ya", "Yb"))
})

test_that("the slopes line names the slopes a class b prior leaves", {
  d <- pts_data()
  fit <- frm(y ~ z + x, data = d, dry_run = "objective")
  line <- function(pr) {
    m <- character()
    withCallingHandlers(pts_rp(fit, pr), message = function(e) {
      m <<- c(m, conditionMessage(e))
      invokeRestart("muffleMessage")
    })
    grep("^  b ", strsplit(paste(m, collapse = "\n"), "\n")[[1]],
         value = TRUE)
  }
  # a class b prior on every slope leaves none flat: no line
  expect_length(line(set_prior("normal(0, 1)", class = "b")), 0L)
  expect_match(line(set_prior("normal(0, 1)", class = "b", coef = "z")),
               "(flat) where your prior sets none", fixed = TRUE)
  # the guard's absent case: no class b prior, the line as it was
  expect_match(line(NULL), "(flat), as brms leaves slopes", fixed = TRUE)
  expect_match(line(set_prior("normal(0, 5)", class = "Intercept")),
               "(flat), as brms leaves slopes", fixed = TRUE)
})

test_that("a uniform on a slope with no bounds is sampled in its support", {
  d <- pts_data()
  fit <- frm(y ~ z, data = d, dry_run = "objective")
  w <- c("no natural lower bound", "no natural upper bound")
  rp <- allow_warnings(suppressMessages(
    pts_rp(fit, set_prior("uniform(0, 0.2)", class = "b", coef = "z"))), w)
  expect_identical(unname(c(rp$ri$lower[["z"]], rp$ri$upper[["z"]])),
                   c(0, 0.2))
})
