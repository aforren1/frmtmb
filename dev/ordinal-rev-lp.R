# Reviewer, lane ordinal: log-density identity against brms 2.23.0's
# compiled Stan program on shapes the worker did not try. Uses the test
# helper's translator (helper-brms.R), but not its expectations, so each
# row prints its numbers instead of stopping at the first failure.
# Seeds: data 20261001 (ord_rev_data). Output: dev/ordinal-rev-log-lp.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
sp <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad"
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(sp, "ordrev-stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
cat("frmtmb from", find.package("frmtmb"), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "tests/testthat/helper-brms.R"), envir = env)
lp_tight <- frmtmb_control(grad_tol = 1e-6, restarts = 3)

lp_row <- function(bform, family, data, fit, joint = FALSE, ...) {
  prior <- env$brms_flat_prior(bform, data = data, family = family, ...)
  code <- brms::make_stancode(bform, data = data, family = family,
                              prior = prior, ...)
  sdat <- env$brms_standata(bform, data = data, family = family,
                            prior = prior, ...)
  rtab <- env$brms_ranef_table(bform, data, family, prior, ...)
  mod <- env$brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  pars <- env$stan_pars_from_fit(fit, sdat, code, rtab)
  upars <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, upars, adjust_transform = FALSE,
                        gradient = FALSE)
  ours <- if (joint) {
    -fit$obj$env$f(fit$obj$env$last.par.best) + attr(pars, "logJ")
  } else as.numeric(logLik(fit))
  grad <- rstan::grad_log_prob(sf, upars, adjust_transform = FALSE)
  if (joint) grad <- grad[env$brms_inner_index(sf, pars)]
  list(lp = lp, ours = ours, const = lp - ours,
       rel = abs(lp - ours) / max(1, abs(ours)), grad = max(abs(grad)),
       pars = pars, code = code)
}

ord_rev_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(letters[1:8], n, TRUE)),
                  h = factor(sample(c("p", "q", "r"), n, TRUE)))
  ug <- rnorm(8, 0, 0.3)[as.integer(d$g)]
  u <- stats::rlogis(n) / exp(0.3 * d$z + ug) + 0.8 * d$x
  d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
  u2 <- stats::rlogis(n) + 0.5 * d$x
  d$y2 <- 1L + (u2 > -1) + (u2 > 0) + (u2 > 1)
  # an unused middle category: 3 of 5 never observed
  d$ym <- d$y
  d$ym[d$ym == 3L] <- 2L
  d
}
d <- ord_rev_data(20261001)
cat("table(ym):", table(d$ym), "\n")

rows <- list(
  list(name = "cumulative, disc ~ x + (1 | g) [joint]",
       f = bf(y ~ x, disc ~ x + (1 | g)), fam = cumulative(),
       bf = brms::bf(y ~ x, disc ~ x + (1 | g)), bfam = brms::cumulative(),
       joint = TRUE),
  list(name = "sratio probit, disc ~ 0 + x + (1 | g) [joint]",
       f = bf(y ~ x, disc ~ 0 + x + (1 | g)), fam = sratio("probit"),
       bf = brms::bf(y ~ x, disc ~ 0 + x + (1 | g)),
       bfam = brms::sratio("probit"), joint = TRUE),
  list(name = "cumulative equidistant, y ~ x + (1 | g) [joint]",
       f = bf(y ~ x + (1 | g)), fam = cumulative(threshold = "equidistant"),
       bf = brms::bf(y ~ x + (1 | g)),
       bfam = brms::cumulative(threshold = "equidistant"), joint = TRUE),
  list(name = "cratio logit, cs(z), disc ~ 0 + x",
       f = bf(y ~ x + cs(z), disc ~ 0 + x), fam = cratio(),
       bf = brms::bf(y ~ x + cs(z), disc ~ 0 + x), bfam = brms::cratio()),
  list(name = "cratio probit equidistant, cs(z), disc ~ 0 + x",
       f = bf(y ~ x + cs(z), disc ~ 0 + x),
       fam = cratio("probit", threshold = "equidistant"),
       bf = brms::bf(y ~ x + cs(z), disc ~ 0 + x),
       bfam = brms::cratio("probit", threshold = "equidistant")),
  list(name = "acat probit, cs(z), disc ~ 0 + x",
       f = bf(y ~ x + cs(z), disc ~ 0 + x), fam = acat("probit"),
       bf = brms::bf(y ~ x + cs(z), disc ~ 0 + x),
       bfam = brms::acat("probit")),
  list(name = "sratio sum_to_zero, cs(z), disc ~ 0 + x",
       f = bf(y ~ x + cs(z), disc ~ 0 + x),
       fam = sratio(threshold = "sum_to_zero"),
       bf = brms::bf(y ~ x + cs(z), disc ~ 0 + x),
       bfam = brms::sratio(threshold = "sum_to_zero")),
  list(name = "cratio sum_to_zero, cs(z)",
       f = bf(y ~ x + cs(z)), fam = cratio(threshold = "sum_to_zero"),
       bf = brms::bf(y ~ x + cs(z)),
       bfam = brms::cratio(threshold = "sum_to_zero")),
  list(name = "acat cloglog sum_to_zero, cs(z)",
       f = bf(y ~ x + cs(z)), fam = acat("cloglog", threshold = "sum_to_zero"),
       bf = brms::bf(y ~ x + cs(z)),
       bfam = brms::acat("cloglog", threshold = "sum_to_zero")),
  list(name = "sratio equidistant, thres(gr = h), disc ~ 0 + z",
       f = bf(y | thres(gr = h) ~ x, disc ~ 0 + z),
       fam = sratio(threshold = "equidistant"),
       bf = brms::bf(y | thres(gr = h) ~ x, disc ~ 0 + z),
       bfam = brms::sratio(threshold = "equidistant")),
  list(name = "acat probit equidistant, thres(gr = h)",
       f = bf(y | thres(gr = h) ~ x),
       fam = acat("probit", threshold = "equidistant"),
       bf = brms::bf(y | thres(gr = h) ~ x),
       bfam = brms::acat("probit", threshold = "equidistant")),
  list(name = "cratio equidistant, thres(gr = h)",
       f = bf(y | thres(gr = h) ~ x), fam = cratio(threshold = "equidistant"),
       bf = brms::bf(y | thres(gr = h) ~ x),
       bfam = brms::cratio(threshold = "equidistant")),
  list(name = "cumulative sum_to_zero, thres(gr = h), disc ~ 0 + z",
       f = bf(y | thres(gr = h) ~ x, disc ~ 0 + z),
       fam = cumulative(threshold = "sum_to_zero"),
       bf = brms::bf(y | thres(gr = h) ~ x, disc ~ 0 + z),
       bfam = brms::cumulative(threshold = "sum_to_zero")),
  list(name = "cumulative equidistant, unused middle category",
       f = bf(ym ~ x), fam = cumulative(threshold = "equidistant"),
       bf = brms::bf(ym ~ x),
       bfam = brms::cumulative(threshold = "equidistant")),
  list(name = "sratio equidistant, unused middle, disc ~ 0 + z",
       f = bf(ym ~ x, disc ~ 0 + z), fam = sratio(threshold = "equidistant"),
       bf = brms::bf(ym ~ x, disc ~ 0 + z),
       bfam = brms::sratio(threshold = "equidistant")),
  list(name = "acat probit equidistant, unused middle",
       f = bf(ym ~ x), fam = acat("probit", threshold = "equidistant"),
       bf = brms::bf(ym ~ x),
       bfam = brms::acat("probit", threshold = "equidistant")),
  list(name = "hurdle_cumulative equidistant, disc ~ 0 + z",
       f = bf(ym ~ x, disc ~ 0 + z),
       fam = hurdle_cumulative(threshold = "equidistant"),
       bf = brms::bf(ym ~ x, disc ~ 0 + z),
       bfam = brms::hurdle_cumulative(threshold = "equidistant")),
  list(name = "mv: cumulative disc ~ 0 + z + sratio sum_to_zero",
       f = bf(y ~ x, disc ~ 0 + z) + bf(y2 ~ x),
       fam = list(cumulative(), sratio(threshold = "sum_to_zero")),
       bf = brms::bf(y ~ x, disc ~ 0 + z) + brms::bf(y2 ~ x) +
         brms::set_rescor(FALSE),
       bfam = list(brms::cumulative(),
                   brms::sratio(threshold = "sum_to_zero"))),
  list(name = "mv: acat probit + cratio (flexible, disc on both)",
       f = bf(y ~ x, disc ~ 0 + z) + bf(y2 ~ x, disc ~ 0 + z),
       fam = list(acat("probit"), cratio()),
       bf = brms::bf(y ~ x, disc ~ 0 + z) + brms::bf(y2 ~ x, disc ~ 0 + z) +
         brms::set_rescor(FALSE),
       bfam = list(brms::acat("probit"), brms::cratio())),
  list(name = "mv: cumulative equidistant (one response only) + sratio",
       f = bf(y ~ x) + bf(y2 ~ x),
       fam = list(cumulative(threshold = "equidistant"), sratio()),
       bf = brms::bf(y ~ x) + brms::bf(y2 ~ x) + brms::set_rescor(FALSE),
       bfam = list(brms::cumulative(threshold = "equidistant"),
                   brms::sratio()))
)

res <- data.frame(row = character(), logLik = numeric(), const = numeric(),
                  rel = numeric(), grad = numeric(), stringsAsFactors = FALSE)
for (r in rows) {
  cat("\n==", r$name, "==\n")
  out <- tryCatch({
    fit <- withCallingHandlers(
      frm(r$f, family = r$fam, data = d, control = lp_tight),
      warning = function(w) {
        cat("frmtmb WARNING:", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")
      })
    chk <- lp_row(r$bf, r$bfam, d, fit, joint = isTRUE(r$joint))
    cat(sprintf("LP const %.6g rel %.3g grad %.3g ours %.10g\n",
                chk$const, chk$rel, chk$grad, chk$ours))
    c(logLik = chk$ours, const = chk$const, rel = chk$rel, grad = chk$grad)
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    c(logLik = NA, const = NA, rel = NA, grad = NA)
  })
  res[nrow(res) + 1L, ] <- list(r$name, out[["logLik"]], out[["const"]],
                                out[["rel"]], out[["grad"]])
}
cat("\n== SUMMARY (data seed 20261001, n = 400) ==\n")
print(format(res, digits = 4), row.names = FALSE)
