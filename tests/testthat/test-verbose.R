# frmtmb_control(verbose =) / frm(verbose =): stage progress and
# timings through message(), and the optimizer's own trace at level 2.

sleep_bf <- function() {
  data(sleepstudy, package = "lme4", envir = environment())
  list(bform = bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
       data = sleepstudy)
}

test_that("a verbose fit reports every stage and a silent one is quiet", {
  s <- sleep_bf()
  msgs <- testthat::capture_messages(
    frm(s$bform, s$data, se = TRUE, verbose = TRUE)
  )
  expect_true(any(grepl("^frmtmb: parse \\[", msgs)))
  expect_true(any(grepl("^frmtmb: frame \\[.*180 obs", msgs)))
  expect_true(any(grepl("^frmtmb: fit: gaussian, ML", msgs)))
  expect_true(any(grepl("^frmtmb: tape \\[.*outer", msgs)))
  expect_true(any(grepl("^frmtmb: optimize \\[.*objective", msgs)))
  expect_true(any(grepl("^frmtmb: sdreport \\[", msgs)))
  expect_true(any(grepl("^frmtmb: done \\[.*max\\|grad\\|.*0 warnings",
                        msgs)))
  # every stage line carries elapsed seconds, so a slow fit localizes
  # itself; only the opening "fit:" line has no timing
  stages <- grep("^frmtmb: fit: ", msgs, invert = TRUE, value = TRUE)
  expect_true(all(grepl("\\[[0-9.]+s\\]", stages)))

  expect_no_message(frm(s$bform, s$data))
  expect_no_message(frm(s$bform, s$data,
                        control = frmtmb_control(verbose = FALSE)))
})

test_that("verbose does not change the fit", {
  s <- sleep_bf()
  f0 <- frm(s$bform, s$data)
  f1 <- suppressMessages(frm(s$bform, s$data, verbose = TRUE))
  # level 2 also writes nlminb's own trace to stdout; keep it out of
  # the test log
  utils::capture.output(
    f2 <- suppressMessages(frm(s$bform, s$data,
                               control = frmtmb_control(verbose = 2)))
  )
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(f0)))
  expect_equal(unlist(fixef(f1)), unlist(fixef(f0)))
})

test_that("control$verbose wins over frm(verbose =)", {
  s <- sleep_bf()
  expect_no_message(frm(s$bform, s$data, verbose = TRUE,
                        control = frmtmb_control(verbose = FALSE)))
  on_by_control <- testthat::capture_messages(
    frm(s$bform, s$data, verbose = FALSE,
        control = frmtmb_control(verbose = TRUE))
  )
  expect_true(any(grepl("^frmtmb: parse", on_by_control)))
  # unset in control, so the frm() shortcut applies
  on_by_frm <- testthat::capture_messages(
    frm(s$bform, s$data, verbose = TRUE,
        control = frmtmb_control(restarts = 0))
  )
  expect_true(any(grepl("^frmtmb: parse", on_by_frm)))
})

test_that("the mode and the autoscale pre-fit are named", {
  s <- sleep_bf()
  dd <- s$data
  dd$Big <- dd$Days * 1e6
  msgs <- testthat::capture_messages(
    frm(bf(Reaction ~ Big + (1 | Subject)) + gaussian(), dd, REML = TRUE,
        control = frmtmb_control(autoscale = TRUE, verbose = TRUE))
  )
  expect_true(any(grepl("fit: gaussian, REML, autoscale", msgs)))
  expect_true(any(grepl("^frmtmb: autoscale pre-fit \\[", msgs)))
  # the pre-fit is one line, not a nested copy of every stage
  expect_length(grep("^frmtmb: tape ", msgs), 1L)

  pmsgs <- testthat::capture_messages(
    frm(s$bform, s$data,
        control = frmtmb_control(profile = TRUE, verbose = TRUE))
  )
  expect_true(any(grepl("fit: gaussian, ML, profile", pmsgs)))
})

test_that("restarts and convergence warnings are reported", {
  s <- sleep_bf()
  msgs <- suppressWarnings(testthat::capture_messages(
    frm(s$bform, s$data,
        control = frmtmb_control(verbose = TRUE, grad_tol = 1e-12,
                                 restarts = 2))
  ))
  expect_length(grep("^frmtmb: restart ", msgs), 2L)
  expect_true(any(grepl("restart 1 \\[.*from max\\|grad\\|", msgs)))
  expect_true(any(grepl("done \\[.*, 1 warning", msgs)))
})

test_that("level 2 sets the optimizer trace but never overrides the user", {
  ctl <- frmtmb_control(verbose = 2)
  expect_equal(vb_trace_ctrl(ctl$optCtrl, "nlminb")$trace, 1L)
  expect_equal(vb_trace_ctrl(ctl$optCtrl, "optim")$REPORT, 1L)
  # an explicit trace is left alone, including trace = 0
  user <- list(iter.max = 1000, trace = 0)
  expect_identical(vb_trace_ctrl(user, "nlminb"), user)
  expect_identical(vb_trace_ctrl(user, "optim"), user)
  # custom optimizers get optCtrl untouched (verbose does not reach them)
  expect_identical(vb_trace_ctrl(ctl$optCtrl, identity), ctl$optCtrl)

  s <- sleep_bf()
  seen <- NULL
  my_opt <- function(par, fn, gr, lower, upper, control) {
    seen <<- control
    r <- stats::nlminb(par, fn, gr, lower = lower, upper = upper)
    list(par = r$par, objective = r$objective,
         convergence = r$convergence)
  }
  suppressMessages(
    frm(s$bform, s$data,
        control = frmtmb_control(optimizer = my_opt, verbose = 2))
  )
  expect_null(seen$trace)
})

test_that("level 1 leaves the optimizer trace off", {
  s <- sleep_bf()
  out <- utils::capture.output(
    invisible(suppressMessages(frm(s$bform, s$data, verbose = TRUE)))
  )
  expect_length(out, 0L)
})

test_that("refit loops stay quiet under a verbose fit", {
  s <- sleep_bf()
  f <- suppressMessages(
    frm(bf(Reaction ~ Days + (1 | Subject)) + gaussian(), s$data,
        control = frmtmb_control(verbose = TRUE))
  )
  expect_no_message(frm_bootstrap(f, nsim = 2, seed = 1))
  expect_no_message(influence(f, groups = "Subject"))
  expect_no_message(suppressWarnings(
    frm_allfit(f, optimizers = list(nlminb = "nlminb"))
  ))
  # a single user-driven refit still reports
  rmsgs <- testthat::capture_messages(refit(f, s$data$Reaction))
  expect_true(any(grepl("^frmtmb: fit: ", rmsgs)))
})

test_that("verbose_level normalizes its input", {
  expect_equal(verbose_level(list()), 0L)
  expect_equal(verbose_level(list(verbose = FALSE)), 0L)
  expect_equal(verbose_level(list(verbose = TRUE)), 1L)
  expect_equal(verbose_level(list(verbose = 2)), 2L)
  expect_equal(verbose_level(list(verbose = -1)), 0L)
  expect_equal(verbose_level(list(verbose = "no")), 0L)
})
