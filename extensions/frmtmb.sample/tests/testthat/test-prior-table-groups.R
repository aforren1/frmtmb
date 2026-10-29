# The group rows of a route = "sample" default_prior() table print
# brms's "" (as "(flat)"), so a class row edited in the table reaches
# every group when the table is passed back to frm_sample(), as in
# brms. Through frmtmb.sample 0.24.0 the group rows of sd, cor and the
# new corme printed the inherited default, which the table then applied
# as an explicit prior over the edit: an edited lkj(4) class row
# sampled lkj(1), KS against Beta(4, 4) p = 1.3e-32
# (dev/reviews/2026-10-28-meme.md, finding 1; brms's rows in
# dev/meme-out/fix-brms-grp.txt).

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

ptg_data <- function(seed = 71) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)),
                  th = factor(rep(c("a", "b"), 30)))
  d$xg <- stats::rnorm(6)[d$g]
  d$zg <- 0.5 * d$xg + stats::rnorm(6, 0, 0.7)[d$g]
  d$sg <- 0.3
  d$y <- 1 + 0.5 * d$x + 0.4 * d$xg + stats::rnorm(6, 0, 0.5)[d$g] +
    stats::rnorm(60)
  d$yo <- cut(d$y, c(-Inf, 0, 1, 2, Inf), labels = FALSE)
  d
}

ptg_form <- function() {
  bf(y ~ x + (1 + x | g) + me(xg, sg, gr = g) + me(zg, sg, gr = g))
}

ptg_rows <- function(tab) {
  tab <- as.data.frame(tab)
  stats::setNames(ifelse(nzchar(tab$prior), tab$prior, "(flat)"),
                  paste(tab$class, tab$coef, tab$group, sep = "|"))
}

ptg_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("the group rows print brms's (flat), as brms prints them", {
  d <- ptg_data()
  got <- ptg_rows(default_prior(ptg_form(), data = d, route = "sample"))
  for (k in c("sd||g", "sd|Intercept|g", "sd|x|g", "cor||g", "corme||g")) {
    expect_identical(got[[k]], "(flat)", label = k)
  }
  expect_identical(got[["sd||"]], "student_t(3, 0, 2.5)")
  expect_identical(got[["cor||"]], "lkj(1)")
  expect_identical(got[["corme||"]], "lkj(1)")
  # grouped thresholds: a level of thres(gr = ) reads "" in brms too
  go <- ptg_rows(default_prior(bf(yo | thres(gr = th) ~ x), data = d,
                               family = cumulative(), route = "sample"))
  expect_identical(go[["Intercept||a"]], "(flat)")
  expect_identical(go[["Intercept||"]], "student_t(3, 0, 2.5)")
  skip_if_not_installed("brms")
  bt <- local({
    attach(list(me = brms::me), name = "ptg_brms_me",
           warn.conflicts = FALSE)
    on.exit(detach("ptg_brms_me"))
    brms::default_prior(brms::bf(y ~ x + (1 + x | g) + me(xg, sg, gr = g) +
                                   me(zg, sg, gr = g)), data = d)
  })
  want <- ptg_rows(bt)
  # every row of brms's table, density by density
  expect_identical(got[names(want)], want)
  bo <- ptg_rows(brms::default_prior(brms::bf(yo | thres(gr = th) ~ x),
                                     data = d,
                                     family = brms::cumulative()))
  expect_identical(go[names(bo)], bo)
})

test_that("an edited class row reaches every group, sd, cor and corme", {
  d <- ptg_data()
  fit <- frm(ptg_form(), data = d, dry_run = "objective")
  tab <- default_prior(ptg_form(), data = d, route = "sample")
  cls <- !nzchar(tab$coef) & !nzchar(tab$group)
  tab$prior[tab$class == "corme" & cls] <- "lkj(4)"
  tab$prior[tab$class == "cor" & cls] <- "lkj(4)"
  tab$prior[tab$class == "sd" & cls] <- "exponential(3)"
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, tab))$ri
  comp <- vapply(ri$entries, `[[`, "", "comp")
  th <- ri$entries[comp == "thetame"]
  expect_length(th, 1L)
  expect_identical(th[[1L]]$dist$eta, 4)
  kinds <- vapply(ri$entries[comp == "theta"], function(e) e$dist$kind, "")
  etas <- unlist(lapply(ri$entries[comp == "theta"], function(e) e$dist$eta))
  # the block's two sds and its correlation
  expect_setequal(kinds, c("exponential", "exponential", "lkj"))
  expect_identical(etas, 4)
})

test_that("the edited class rows are what the prior draws come from", {
  skip_sampler()
  d <- ptg_data()
  tab <- default_prior(ptg_form(), data = d, route = "sample")
  cls <- !nzchar(tab$coef) & !nzchar(tab$group)
  tab$prior[tab$class == "corme" & cls] <- "lkj(4)"
  tab$prior[tab$class == "cor" & cls] <- "lkj(4)"
  tab$prior[tab$class == "sd" & cls] <- "exponential(3)"
  ds <- ptg_quiet(frm_sample(ptg_form(), data = d, family = gaussian(),
                             prior = tab, sample_prior = "yes", chains = 2,
                             iter = 1000, seed = 4, refresh = 0))
  pd <- prior_draws(ds)
  ks <- function(v, cdf, ...) {
    suppressWarnings(stats::ks.test(v, cdf, ...))$p.value
  }
  # one correlation of an LKJ(4) 2 x 2 matrix is 2 Beta(4, 4) - 1
  expect_gt(ks((pd$Lme_1 + 1) / 2, "pbeta", 4, 4), 0.01)
  expect_gt(ks((pd$cor_g + 1) / 2, "pbeta", 4, 4), 0.01)
  expect_gt(ks(pd$sd_g, "pexp", 3), 0.01)
  # and not the lkj(1) the group rows used to impose
  expect_lt(ks((pd$Lme_1 + 1) / 2, "pbeta", 1, 1), 0.01)
})

test_that("without me(), the sd and cor group rows round-trip too", {
  # the same trap on a model 0.24.0 could already sample: the edited
  # lkj(4) and exponential(3) class rows sampled lkj(1) and the
  # student_t of the printed group row (dev/review-meme-out/grpinh2.txt)
  d <- ptg_data()
  f <- bf(y ~ x + (1 + x | g))
  fit <- frm(f, data = d, dry_run = "objective")
  tab <- default_prior(f, data = d, route = "sample")
  got <- ptg_rows(tab)
  expect_identical(unname(got[c("sd||g", "cor||g")]), c("(flat)", "(flat)"))
  cls <- !nzchar(tab$coef) & !nzchar(tab$group)
  tab$prior[tab$class == "cor" & cls] <- "lkj(4)"
  tab$prior[tab$class == "sd" & cls] <- "exponential(3)"
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, tab))$ri
  th <- Filter(function(e) identical(e$comp, "theta"), ri$entries)
  expect_setequal(vapply(th, function(e) e$dist$kind, ""),
                  c("exponential", "exponential", "lkj"))
  expect_identical(unlist(lapply(th, function(e) e$dist$eta)), 4)
})
