# Short check of fix round 2, item 3: F2's rule through update(), sample_prior and
# at ordinary noise. Which init did each call use, read off the stanfit's own start?
# Does a good mode now start at random, and what does that cost? Seed 6601.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-f2.R f2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(6601)
n <- 15
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + 1e-3 * e0)
f0 <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
cat(sprintf("frm() sigma %.3g\n", as.numeric(sigma(f0))))
sigstart <- function(fs) {
  ii <- rstan::get_inits(fs$stanfit)
  v <- vapply(ii, function(z) as.numeric(unlist(z))[2L], 0)
  paste(signif(v, 3), collapse = ",")
}
run <- function(tag, expr) {
  msg <- character()
  fs <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      cat(sprintf("F2 %-34s ERROR %s\n", tag, substr(conditionMessage(e), 1, 120)))
      NULL
    }),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    }, warning = function(w) invokeRestart("muffleWarning"))
  if (is.null(fs)) return(invisible(NULL))
  said <- any(grepl("chains start at random", msg, fixed = TRUE))
  dd <- posterior::as_draws_df(fs)
  st <- tryCatch(sigstart(fs), error = function(e) "n/a")
  cat(sprintf("F2 %-34s says random %-5s | start betad (log sigma): %s | post-warmup sigma median %.3g\n",
              tag, said, st, median(dd$sigma)))
  invisible(fs)
}
a <- list(chains = 2, iter = 120, warmup = 60, seed = 5, refresh = 0, cores = 1)
fs <- run("frm_sample(fit)", do.call(frm_sample, c(list(f0), a)))
run("frm_sample(fit, last.par.best)", do.call(frm_sample, c(list(f0, init = "last.par.best"), a)))
run("update(draws)", update(fs, refresh = 0))
run("update(draws, iter = 100)", update(fs, iter = 100, warmup = 50, refresh = 0))
run("update(draws, new prior)", update(fs, prior = set_prior("normal(0, 1)", class = "sdgp"), refresh = 0))
run("update(draws, newdata)", update(fs, newdata = d[1:12, ], refresh = 0))
run("frm_sample(fit, sample_prior = yes)", do.call(frm_sample, c(list(f0, sample_prior = "yes"), a)))
run("frm_sample(fit, sample_prior = only)", do.call(frm_sample, c(list(f0, sample_prior = "only"), a)))
run("frm_sample(fit, reparameterize = F)", do.call(frm_sample, c(list(f0, reparameterize = FALSE), a)))
run("frm_sample(formula)", do.call(frm_sample, c(list(bf(y ~ gp(x)), data = d, family = gaussian()), a)))
d2 <- transform(d, y = sin(2 * x) + 0.1 * e0)
f1 <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d2))
run("ordinary noise frm_sample(fit)", do.call(frm_sample, c(list(f1), a)))
d3 <- transform(d, z = sin(x))
f3 <- suppressWarnings(frm(bf(y ~ z), family = gaussian(), data = d3))
run("no gp frm_sample(fit)", do.call(frm_sample, c(list(f3), a)))
f4 <- suppressWarnings(frm(bf(y ~ gp(x, k = 10)), family = gaussian(), data = d))
run("hsgp frm_sample(fit)", do.call(frm_sample, c(list(f4), a)))
f5 <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = gaussian(),
                           data = transform(d, g = factor(rep(1:3, length.out = n)))))
run("(1|g) frm_sample(fit)", do.call(frm_sample, c(list(f5), a)))

# the cost at ordinary noise: mode against random on a good fit
set.seed(6611)
nn <- 60
xx <- round(runif(nn, 0, 6), 3)
dd <- data.frame(x = xx, y = sin(2 * xx) + 0.3 * rnorm(nn))
ff <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = dd))
for (ini in c("random", "last.par.best")) {
  t0 <- proc.time()[["elapsed"]]
  fs2 <- suppressMessages(suppressWarnings(frm_sample(ff, init = ini, chains = 4, iter = 1000, warmup = 500,
                                                      seed = 7, refresh = 0, cores = 1)))
  w <- proc.time()[["elapsed"]] - t0
  dv <- posterior::as_draws_df(fs2)
  m <- function(v) posterior::extract_variable_matrix(dv, v)
  sp <- rstan::get_sampler_params(fs2$stanfit, inc_warmup = FALSE)
  cat(sprintf("F2COST n 60 noise 0.3 init %-14s wall %.0f s | Rhat sigma %.3f sdgp %.3f lscale %.3f | ESS %.0f %.0f %.0f | divergent %d\n",
              ini, w, posterior::rhat(m("sigma")), posterior::rhat(m("sdgp_gpx")),
              posterior::rhat(m("lscale_gpx")), posterior::ess_bulk(m("sigma")),
              posterior::ess_bulk(m("sdgp_gpx")), posterior::ess_bulk(m("lscale_gpx")),
              sum(vapply(sp, function(s) sum(s[, "divergent__"]), 0))))
}
