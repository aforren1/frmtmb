# brms's class "sds", the standard deviation of a smooth's penalized part:
# its spelling, its reach and its rows on default_prior(), as brms 2.23.0
# has them (dev/refuse-brms-sds.R, dev/refuse-table-compare.R).

sds_data <- function(seed = 43, n = 120) {
  set.seed(seed)
  d <- data.frame(x = runif(n), z = runif(n), w = rnorm(n),
                  f = factor(rep(c("a", "b"), length.out = n)),
                  g = factor(rep(1:10, length.out = n)))
  d$y <- sin(2 * pi * d$x) + d$z + rnorm(n, 0, 0.3)
  d$y2 <- d$x + rnorm(n, 0, 0.3)
  d
}

# what a prior resolves to on a design, without a fit
sds_entries <- function(form, pr, dd) {
  fr <- frm(form, data = dd, dry_run = "frame")
  sp <- frm(form, data = dd, dry_run = "spec")
  frmtmb:::resolve_priorlist(list(spec = sp, frame = fr), pr)$entries
}

test_that("class sds is one of set_prior()'s classes, on both spellings", {
  p <- set_prior("student_t(3, 0, 1)", class = "sds")
  expect_identical(unclass(p)[[1L]]$class, "sds")
  expect_null(frmtmb:::brms_prior_class_refusal("sds"))
  expect_identical(frmtmb:::brms_prior_route("sds", "normal(0, 1)")$class,
                   "sds")
  skip_if_not_installed("brms")
  pl <- frmtmb:::as_priorlist(brms::prior(student_t(3, 0, 1), class = "sds"))
  expect_identical(unclass(pl)[[1L]]$class, "sds")
  # a term is named by coef; a smooth is not a grouping term
  expect_error(set_prior("normal(0, 1)", class = "sds", group = "s(x)"),
               "has no group")
})

test_that("default_prior() lists a smooth's rows as brms does", {
  skip_if_not_installed("brms")
  dd <- sds_data()
  forms <- list(
    list(bf(y ~ s(x)), brms::bf(y ~ s(x))),
    list(bf(y ~ s(x) + s(z, k = 5)), brms::bf(y ~ s(x) + s(z, k = 5))),
    list(bf(y ~ s(z, bs = "cr", k = 5)), brms::bf(y ~ s(z, bs = "cr", k = 5))),
    list(bf(y ~ t2(x, z)), brms::bf(y ~ t2(x, z))),
    list(bf(y ~ f + s(x, by = f)), brms::bf(y ~ f + s(x, by = f))),
    list(bf(y ~ s(x) + (1 | g)), brms::bf(y ~ s(x) + (1 | g))),
    list(bf(y ~ s(x), sigma ~ s(z)), brms::bf(y ~ s(x), sigma ~ s(z))),
    list(bf(y ~ a * exp(b * x), a ~ s(z), b ~ 1, nl = TRUE),
         brms::bf(y ~ a * exp(b * x), a ~ s(z), b ~ 1, nl = TRUE)),
    list(mvbf(bf(y ~ s(x)), bf(y2 ~ s(z)), rescor = FALSE),
         brms::bf(y ~ s(x)) + brms::bf(y2 ~ s(z)) + brms::set_rescor(FALSE))
  )
  cols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
  rows <- function(tab) {
    tab <- as.data.frame(tab)
    tab <- tab[tab$class %in% c("sd", "sds"), , drop = FALSE]
    out <- as.data.frame(lapply(tab[, cols], as.character))
    # brms's lb on the sd class row is filed; only sds's is read here
    out$lb <- ifelse(tab$class == "sds" & !is.na(tab$lb) &
                       as.character(tab$lb) != "", as.character(tab$lb), "")
    rownames(out) <- NULL
    out
  }
  for (fm in forms) {
    ours <- rows(default_prior(fm[[1L]], data = dd))
    theirs <- rows(brms::default_prior(fm[[2L]], data = dd))
    expect_equal(ours, theirs)
    # the guard is not empty: every form has a class row and a term row
    expect_true(any(ours$class == "sds" & nzchar(ours$coef)))
    expect_true(any(ours$class == "sds" & ours$lb == "0"))
  }
})

test_that("a model with no smooth has no sds row", {
  dd <- sds_data()
  tab <- as.data.frame(default_prior(bf(y ~ x + (1 | g)), data = dd))
  expect_false("sds" %in% tab$class)
  expect_true("sd" %in% tab$class)
})

test_that("class sds reaches every smooth of its prefix, and nothing else", {
  dd <- sds_data()
  pr <- set_prior("normal(0, 1)", class = "sds")
  # one per smooth; the random effect beside it keeps its slot
  e <- sds_entries(bf(y ~ s(x) + s(z) + (1 | g)), pr, dd)
  expect_length(e, 2L)
  # every penalty of a t2() and every level of a by-smooth
  expect_length(sds_entries(bf(y ~ t2(x, z)), pr, dd), 3L)
  expect_length(sds_entries(bf(y ~ f + s(x, by = f)), pr, dd), 2L)
  # no dpar is the location only, as in brms
  expect_length(sds_entries(bf(y ~ s(x), sigma ~ s(z)), pr, dd), 1L)
  expect_length(sds_entries(bf(y ~ s(x), sigma ~ s(z)),
                            set_prior("normal(0, 1)", class = "sds",
                                      dpar = "sigma"), dd), 1L)
  # the density is on the standard deviation itself
  expect_identical(e[[1L]]$scale, "sd")
})

test_that("coef names one term as written, and only that", {
  dd <- sds_data()
  e <- sds_entries(bf(y ~ s(x) + s(z, k = 5)),
                   set_prior("normal(0, 1)", class = "sds",
                             coef = "s(z, k = 5)"),
                   dd)
  expect_length(e, 1L)
  e2 <- sds_entries(bf(y ~ f + s(x, by = f)),
                    set_prior("normal(0, 1)", class = "sds",
                              coef = "s(x, by = f)"), dd)
  expect_length(e2, 2L)
  # brms refuses each of these: "s(z,k=5)" and a by-level are not rows
  expect_error(sds_entries(bf(y ~ s(z, k = 5)),
                           set_prior("normal(0, 1)", class = "sds",
                                     coef = "s(z,k=5)"), dd),
               "No smooth standard deviations match", fixed = TRUE)
  expect_error(sds_entries(bf(y ~ f + s(x, by = f)),
                           set_prior("normal(0, 1)", class = "sds",
                                     coef = "s(x):fa"), dd),
               "s(x, by = f)", fixed = TRUE)
  expect_error(sds_entries(bf(y ~ x),
                           set_prior("normal(0, 1)", class = "sds"), dd),
               "no smooth terms", fixed = TRUE)
})

test_that("class sd reaches no smooth, as in brms", {
  dd <- sds_data()
  pr <- set_prior("normal(0, 1)", class = "sd")
  # beside a random effect it is the random effect's alone
  e <- sds_entries(bf(y ~ s(x) + (1 | g)), pr, dd)
  expect_length(e, 1L)
  # with no random effect, and by the smooth's old group, it names nothing
  expect_error(sds_entries(bf(y ~ s(x)), pr, dd), "class \"sds\"",
               fixed = TRUE)
  expect_error(sds_entries(bf(y ~ s(x) + (1 | g)),
                           set_prior("normal(0, 1)", class = "sd",
                                     group = "s(x)"), dd),
               "class \"sds\"", fixed = TRUE)
})

test_that("the prior is the MAP's and the printed table says so", {
  dd <- sds_data()
  pr <- set_prior("student_t(3, 0, 1)", class = "sds")
  fit <- frm(bf(y ~ s(x)), data = dd, prior = pr)
  plain <- frm(bf(y ~ s(x)), data = dd)
  # the penalty is on sds: the objectives at one internal point differ by
  # the density and the log map's Jacobian, whose sum names no constant
  p <- plain$opt$par
  sd <- exp(p[grepl("^theta", names(p))])
  dens <- stats::dt(sd, 3, log = TRUE) + log(sd)
  expect_lt(abs((plain$obj$fn(p) - fit$obj$fn(p)) - sum(dens)) /
              max(1, abs(sum(dens))), 1e-6)
  ps <- as.data.frame(prior_summary(fit))
  r <- ps[ps$class == "sds", , drop = FALSE]
  expect_true(any(r$prior == "student_t(3, 0, 1)" & !nzchar(r$coef)))
  tab <- validate_prior(pr, bf(y ~ s(x)), data = dd)
  expect_identical(tab$source[tab$class == "sds" & tab$coef == "s(x)"],
                   "(vectorized)")
  # the table goes back in as a prior
  back <- frm(bf(y ~ s(x)), data = dd, prior = tab)
  expect_lt(abs(fit$opt$objective - back$opt$objective) /
              max(1, abs(fit$opt$objective)), 1e-8)
})
