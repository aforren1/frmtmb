# frm_allfit() starts every optimizer where lme4's allFit() does, and
# its table says when a success code sits at a worse optimum. Lane
# nanse (defect 9 of lane vigport), dev/nanse-findings.md.

allfit_nl_fit <- function() {
  set.seed(3)
  d <- data.frame(x = runif(80, 0, 10))
  d$y <- 5 * (1 - exp(-d$x / 3)) + rnorm(80, 0, 0.2)
  frm(bf(y ~ a * (1 - exp(-x / b)), a ~ 1, b ~ 1, nl = TRUE), data = d,
      start = list(beta = c(4, 2)))
}

test_that("a nonlinear fit's refits start at its estimates", {
  # the refits used start = NULL: from zero, b = 0 divides by zero, and
  # 3 of 4 refits returned NULL while bobyqa reported code 0 at a
  # log-likelihood 398 below the fit it was checking
  fit <- allfit_nl_fit()
  af <- frm_allfit(fit)
  expect_true(all(!vapply(af$fits, is.null, TRUE)))
  ll <- vapply(af$fits, function(f) as.numeric(logLik(f)), 0)
  ll0 <- as.numeric(logLik(fit))
  # agreement measured against the fit's own tolerance on log likelihood
  expect_true(all(abs(ll - ll0) < fit$control$grad_tol))
  out <- capture.output(print(af))
  expect_false(any(grepl("converged elsewhere|failed", out)))
})

test_that("start_from_mle = FALSE starts where the original fit did", {
  fit <- allfit_nl_fit()
  af <- frm_allfit(fit, optimizers = list(nlminb = "nlminb"),
                   start_from_mle = FALSE)
  expect_false(is.null(af$fits$nlminb))
  expect_equal(as.numeric(logLik(af$fits$nlminb)),
               as.numeric(logLik(fit)))
  # the user's start is what the refit was given
  expect_identical(af$fits$nlminb$start, fit$start)
  expect_error(frm_allfit(fit, start_from_mle = NA), "start_from_mle")
})

test_that("a success code at a worse optimum is not shown as agreement", {
  fit <- allfit_nl_fit()
  # an optimizer that reports success one unit away from where it began
  stay_off <- function(par, fn, gr, lower, upper, control) {
    p <- par + 0.5
    list(par = p, objective = fn(p), convergence = 0L, message = "ok")
  }
  af <- frm_allfit(fit, optimizers = list(nlminb = "nlminb",
                                          stay_off = stay_off))
  expect_identical(af$fits$stay_off$opt$convergence, 0L)
  out <- capture.output(print(af))
  row <- grep("^ *stay_off", out, value = TRUE)
  expect_match(row, "converged elsewhere", fixed = TRUE)
  expect_match(grep("^ *nlminb", out, value = TRUE), "agrees", fixed = TRUE)
  expect_true(any(grepl("stopped at a different point", out, fixed = TRUE)))
})

test_that("a refit that errors keeps its message", {
  fit <- allfit_nl_fit()
  broken <- function(par, fn, gr, lower, upper, control) stop("no dice")
  af <- frm_allfit(fit, optimizers = list(broken = broken))
  expect_null(af$fits$broken)
  expect_match(af$errors[["broken"]], "no dice", fixed = TRUE)
  expect_match(paste(capture.output(print(af)), collapse = "\n"),
               "failed:", fixed = TRUE)
})
