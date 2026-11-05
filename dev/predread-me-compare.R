# me() on new data, frmtmb against brms 2.23.0 at the same parameter values.
# Part A: frm() and a Fixed_param brmsfit at frm()'s estimate (data seed 41,
# 4000 draws a side). Part B: frmtmb.sample's own draws put through brms:
# each of 200 frm_sample() draws is one chain of a Fixed_param brmsfit, so the
# two packages read the same draws (sampler seed 5).
#   Rscript dev/predread-me-compare.R <lib> > dev/predread-log/me-compare-<tag>.txt
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r16"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r16",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = E)
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
set.seed(41)
n <- 120
g <- rep(1:12, each = 10)
d <- data.frame(g = factor(g), z = rnorm(n))
lat <- rnorm(12)[g]
d$xg <- lat + rnorm(12, 0, 0.3)[g]
d$sg <- 0.3
d$x <- rnorm(n)
d$sx <- 0.4
d$w <- 0.6 * d$x + rnorm(n)
d$sw <- 0.25
d$y <- 1 + 0.8 * d$x + 0.5 * d$w + 0.6 * lat + 0.3 * d$z + rnorm(n, 0, 0.5)
nd <- data.frame(x = c(-1, 0, 1), sx = c(0.2, 0.4, 0.8), w = 0.5, sw = 0.25,
                 z = c(0, 1, -1), xg = c(0.5, 0.5, -1), sg = 0.3,
                 g = factor(c("1", "1", "2"), levels = levels(d$g)))
zcol <- function(a, b, what) {
  za <- (colMeans(a) - colMeans(b)) /
    sqrt(apply(a, 2, var) / nrow(a) + apply(b, 2, var) / nrow(b))
  rv <- apply(a, 2, sd) / apply(b, 2, sd)
  cat(sprintf("  %-24s mean z: %s; sd ratio: %s\n", what,
              paste(sprintf("%.2f", za), collapse = " "),
              paste(sprintf("%.3f", rv), collapse = " ")))
}
models <- list(
  "me(x, sx) + z" = list(f = bf(y ~ me(x, sx) + z),
                         b = brms::bf(y ~ me(x, sx) + z)),
  "me(xg, sg, gr = g) + z" = list(f = bf(y ~ me(xg, sg, gr = g) + z),
                                  b = brms::bf(y ~ me(xg, sg, gr = g) + z)),
  "me(x, sx) + me(w, sw)" = list(f = bf(y ~ me(x, sx) + me(w, sw)),
                                 b = brms::bf(y ~ me(x, sx) + me(w, sw))),
  "me(x, sx) + me(x, sx):z" = list(f = bf(y ~ me(x, sx) + me(x, sx):z),
                                   b = brms::bf(y ~ me(x, sx) + me(x, sx):z)))
cat("\n== A: frm() estimate, Fixed_param brmsfit, 4000 draws a side\n")
for (nm in names(models)) {
  cat(nm, "\n")
  m <- models[[nm]]
  fit <- quiet(frm(m$f, data = d, family = gaussian()))
  bf1 <- E$brms_fixed_fit(m$b, stats::gaussian(), d, fit, ndraws = 4000)
  bp <- quiet(brms::posterior_predict(bf1, newdata = nd))
  fp <- quiet(predict(fit, newdata = nd, summary = FALSE, ndraws = 4000,
                      propagate_error = FALSE))
  zcol(fp, bp, "predict")
  be <- quiet(brms::posterior_epred(bf1, newdata = nd))
  ft <- quiet(fitted(fit, newdata = nd))
  cat(sprintf("  fitted Estimate - brms mean over draws / MCSE: %s\n",
              paste(sprintf("%.2f", (ft[, "Estimate"] - colMeans(be)) /
                              (apply(be, 2, sd) / sqrt(nrow(be)))),
                    collapse = " ")))
  cat(sprintf("  fitted Est.Error %s; coefficient part alone would be %s; ",
              paste(sprintf("%.4f", ft[, "Est.Error"]), collapse = " "),
              paste(sprintf("%.4f", quiet(frm_linpred(
                fit, newdata = nd, dpar = "mu", se.fit = TRUE))$se.fit),
                collapse = " ")))
  cat(sprintf("brms sd at the estimate %s\n",
              paste(sprintf("%.4f", apply(be, 2, sd)), collapse = " ")))
}
cat("\n== B: frm_sample() draws, the same 200 draws in both packages\n")
for (nm in names(models)) {
  cat(nm, "\n")
  m <- models[[nm]]
  ds <- quiet(frm_sample(m$f, data = d, family = gaussian(),
                         chains = 2, iter = 500, seed = 5, refresh = 0))
  draws <- round(seq(1, posterior::ndraws(ds), length.out = 200))
  prior <- E$brms_flat_prior(m$b, data = d, family = stats::gaussian())
  code <- brms::make_stancode(m$b, data = d, family = stats::gaussian(),
                              prior = prior)
  sdat <- E$brms_standata(m$b, data = d, family = stats::gaussian(),
                          prior = prior)
  scaffold <- quiet(brms::brm(m$b, data = d, family = stats::gaussian(),
                              prior = prior, empty = TRUE))
  idx <- frmtmb.sample:::draws_index(ds)
  inits <- lapply(draws, function(k) {
    sh <- frmtmb.sample:::draws_fit_at(ds, k, idx)
    E$stan_pars_from_fit(sh, sdat, code, scaffold$ranef)
  })
  mod <- E$brms_stan_model(code)
  sf <- quiet(rstan::sampling(mod, data = sdat, chains = length(inits),
                              iter = 40, warmup = 0, algorithm = "Fixed_param",
                              init = inits, refresh = 0, seed = 1, cores = 1))
  scaffold$fit <- sf
  bfit <- brms::rename_pars(scaffold)
  be <- quiet(brms::posterior_epred(bfit, newdata = nd))
  bp <- quiet(brms::posterior_predict(bfit, newdata = nd))
  bl <- quiet(brms::posterior_linpred(bfit, newdata = nd))
  set.seed(9)
  rep_draw <- function(f) do.call(rbind, lapply(seq_len(40), function(i) {
    quiet(f(ds, newdata = nd, draw_ids = draws))
  }))
  fe <- rep_draw(posterior_epred)
  fp <- rep_draw(posterior_predict)
  fl <- rep_draw(posterior_linpred)
  zcol(fe, be, "posterior_epred")
  zcol(fp, bp, "posterior_predict")
  zcol(fl, bl, "posterior_linpred")
  if (grepl("gr = g", nm)) {
    cat(sprintf("  rows 1 and 2 (one level): frmtmb identical %s; brms %s\n",
                identical(fe[, 1] - fe[, 2], rep(fe[1, 1] - fe[1, 2],
                                                 nrow(fe))) ||
                  sd(fe[, 1] - fe[, 2]) < 1e-1,
                sd(be[, 1] - be[, 2]) < 1e-1))
    cat(sprintf("  sd(row 1 - row 2): frmtmb %.4f brms %.4f\n",
                sd(fe[, 1] - fe[, 2]), sd(be[, 1] - be[, 2])))
    cat(sprintf("  sd(row 1 - row 3): frmtmb %.4f brms %.4f\n",
                sd(fe[, 1] - fe[, 3]), sd(be[, 1] - be[, 3])))
  }
}
