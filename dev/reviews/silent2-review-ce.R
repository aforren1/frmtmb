# Review of lane silent2, item 5362: conditional_effects(method =
# "posterior_predict") on laplace draws across families; what used to run
# and now refuses or errors, and the cost of the extra dpar reads.
#   Rscript dev/reviews/silent2-review-ce.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-silent2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
lap_pair <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, setdiff(colnames(M), inner), drop = FALSE]
  list(full = full, lap = lap)
}
try_s <- function(e) {
  ws <- character()
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(withCallingHandlers({force(e); "ran"},
                  warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") },
                  message = function(m) invokeRestart("muffleMessage")),
                error = function(err) paste("ERROR:", substr(gsub("\n", " ", conditionMessage(err)), 1, 150)))
  el <- proc.time()[["elapsed"]] - t0
  paste0(r, sprintf(" (%.2fs)", el),
         if (length(ws)) paste0(" [", length(ws), " warn: ", substr(ws[1], 1, 90), "]"))
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
set.seed(5)
n <- 120
d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 20)))
d$y <- 0.3 * d$x + rnorm(6)[d$h] + rnorm(n)
d$y2 <- d$x + rnorm(n)
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x)) * rbinom(n, 1, 0.7)
d$pos <- rgamma(n, 2, 2 / exp(0.2 * d$x)) * rbinom(n, 1, 0.8)
d$o <- cut(d$y, c(-Inf, -1, 0, 1, Inf), labels = FALSE)
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$yb <- rbinom(n, 1, plogis(d$x))
d$ynl <- 2 * exp(-0.5 * d$x) + rnorm(n, 0, 0.2)
models <- list(
  gauss_sigma_x = list(bf(y ~ x + (1 | h), sigma ~ x), gaussian()),
  student = list(bf(y ~ x + (1 | h)), student()),
  zip_zi_x = list(bf(cnt ~ x + (1 | h), zi ~ x), zero_inflated_poisson()),
  hurdle_gamma = list(bf(pos ~ x + (1 | h), hu ~ x), hurdle_gamma()),
  mix_gauss = list(bf(y ~ x + (1 | h)), mixture(gaussian(), gaussian())),
  categorical = list(bf(cat ~ x + (1 | h)), categorical()),
  cumulative = list(bf(o ~ x + (1 | h)), cumulative()),
  cumulative_disc = list(bf(o ~ x + (1 | h), disc ~ x), cumulative()),
  bernoulli = list(bf(yb ~ x + (1 | h)), bernoulli()),
  nonlinear = list(bf(ynl ~ a * exp(b * x), a ~ 1 + (1 | h), b ~ 1, nl = TRUE), gaussian()),
  mv = list(bf(y ~ x + (1 | h)) + bf(y2 ~ x) + set_rescor(FALSE), gaussian()),
  gauss_no_re = list(bf(y ~ x, sigma ~ x), gaussian()),
  gauss_sigma_re = list(bf(y ~ x, sigma ~ (1 | h)), gaussian())
)
for (nm in names(models)) {
  m <- models[[nm]]
  fit <- tryCatch(q(m[[1]], family = m[[2]], data = d,
                    prior = if (nm == "nonlinear") c(set_prior("normal(2, 1)", nlpar = "a"),
                                                       set_prior("normal(0, 1)", nlpar = "b")) else NULL),
                  error = function(e) e)
  if (inherits(fit, "error")) { cat(sprintf("%-18s FIT ERROR %s\n", nm, conditionMessage(fit))); next }
  p <- tryCatch(lap_pair(fit), error = function(e) NULL)
  if (is.null(p)) { cat(sprintf("%-18s lap_pair helper does not apply\n", nm)); next }
  resp <- if (nm == "mv") "y" else NULL
  set.seed(4)
  cat(sprintf("%-18s lap re_formula=NA   %s\n", nm,
              try_s(conditional_effects(p$lap, effects = "x", re_formula = NA, resp = resp,
                                        method = "posterior_predict"))))
  set.seed(4)
  cat(sprintf("%-18s lap re_formula=NULL %s\n", nm,
              try_s(conditional_effects(p$lap, effects = "x", re_formula = NULL, resp = resp,
                                        method = "posterior_predict"))))
  set.seed(4)
  cat(sprintf("%-18s lap posterior_epred %s\n", nm,
              try_s(conditional_effects(p$lap, effects = "x", re_formula = NA, resp = resp))))
}
# cost: 200 laplace draws, gaussian with sigma ~ x, re_formula = NA
fit <- q(bf(y ~ x + (1 | h), sigma ~ x), data = d)
p <- lap_pair(fit, n = 200L)
t <- system.time(conditional_effects(p$lap, effects = "x", re_formula = NA,
                                      method = "posterior_predict"))
cat("COST 200 draws sigma ~ x, re_formula = NA: elapsed", t[["elapsed"]], "\n")
cat("DONE\n")
