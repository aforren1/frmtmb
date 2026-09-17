## Recheck (rounds 1 and 1b), priority 3: the trials() refusals for
## binomial, beta_binomial, zero_inflated_binomial, their mixtures and
## multinomial, on the response shapes the field writes, in three arms:
## brms (make_standata / get_prior), base and lane (frm, and the
## post-fit paths). Seed 20260916, n = 120.
## Usage: Rscript dev/famlink-rev2-trials.R <lane|base|brms>
ARM <- commandArgs(trailingOnly = TRUE)[1]
if (identical(ARM, "brms")) {
  .libPaths(c("C:/Users/adf44/source/r/pinlib",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(brms))
} else {
  source("dev/famlink-rev-common.R")
}
set.seed(20260916)
n <- 120
d <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)))
d$n <- sample(4:12, n, TRUE)
d$s <- rbinom(n, d$n, plogis(0.2 + 0.5 * d$x))
d$f <- d$n - d$s
d$y01 <- rbinom(n, 1, plogis(0.3 * d$x))
d$prop <- d$s / d$n
d$w <- runif(n, 0.5, 2)
d$one <- 1L
d$y2 <- rbinom(n, d$n, 0.4)
P <- t(sapply(seq_len(n), function(i) rmultinom(1, d$n[i], c(0.2, 0.3, 0.5))))
d$Y <- P
d$Yone <- t(sapply(seq_len(n), function(i) rmultinom(1, 1, c(0.2, 0.3, 0.5))))
F <- if (identical(ARM, "brms")) brms::bf else frmtmb::bf
mvb <- if (identical(ARM, "brms")) brms::mvbf else frmtmb::mvbf
MIX <- if (identical(ARM, "brms")) brms::mixture else frmtmb::mixture
MULTI <- function(K) if (identical(ARM, "brms")) brms::multinomial() else frmtmb::multinomial(K)
BIN <- function(...) if (identical(ARM, "brms")) brms::brmsfamily("binomial", ...) else stats::binomial(...)
BB <- function() if (identical(ARM, "brms")) brms::beta_binomial() else frmtmb::beta_binomial()
ZIB <- function() if (identical(ARM, "brms")) brms::zero_inflated_binomial() else frmtmb::zero_inflated_binomial()
cases <- list(
  trials_col = list(s | trials(n) ~ x, quote(BIN())),
  trials_const1_y01 = list(y01 | trials(1) ~ x, quote(BIN())),
  no_trials_y01 = list(y01 ~ x, quote(BIN())),
  no_trials_counts = list(s ~ x, quote(BIN())),
  cbind_sf = list(cbind(s, f) ~ x, quote(BIN())),
  cbind_sf_re = list(cbind(s, f) ~ x + (1 | g), quote(BIN())),
  cbind_betabin = list(cbind(s, f) ~ x, quote(BB())),
  weights_y01_no_trials = list(y01 | weights(w) ~ x, quote(BIN())),
  prop_weights_no_trials = list(prop | weights(n) ~ x, quote(BIN())),
  trials_and_weights = list(s | trials(n) + weights(w) ~ x, quote(BIN())),
  betabin_no_trials = list(y01 ~ x, quote(BB())),
  zib_no_trials = list(y01 ~ x, quote(ZIB())),
  zib_trials = list(s | trials(n) ~ x, quote(ZIB())),
  family_string_no_trials = list(y01 ~ x, "binomial"),
  family_string_trials = list(s | trials(n) ~ x, "binomial"),
  mv_both_trials = list(quote(mvb(F(s | trials(n) ~ x), F(y2 | trials(n) ~ x)) + BIN()), NULL),
  mv_one_missing = list(quote(mvb(F(s | trials(n) ~ x), F(y01 ~ x)) + BIN()), NULL),
  mv_mixed_families = list(quote(mvb(F(s | trials(n) ~ x) + BIN(), F(x ~ 1) + gaussian())), NULL),
  mix_binom_trials = list(s | trials(n) ~ 1, quote(MIX(BIN(), BIN()))),
  mix_binom_no_trials = list(y01 ~ 1, quote(MIX(BIN(), BIN()))),
  mix_binom_poisson_trials = list(s | trials(n) ~ 1, quote(MIX(BIN(), poisson()))),
  mix_bern_bern = list(y01 ~ 1, quote(MIX(bernoulli(), bernoulli()))),
  multinom_trials = list(Y | trials(n) ~ x, quote(MULTI(3))),
  multinom_no_trials = list(Y ~ x, quote(MULTI(3))),
  multinom_trials_wrong = list(Y | trials(one) ~ x, quote(MULTI(3))),
  multinom_onehot_trials1 = list(Yone | trials(1) ~ x, quote(MULTI(3))),
  multinom_onehot_no_trials = list(Yone ~ x, quote(MULTI(3)))
)
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fo <- cs[[1]]
  if (is.call(fo) && !inherits(fo, "formula")) fo <- tryCatch(eval(fo), error = function(e) e)
  fam <- if (is.null(cs[[2]])) NULL else if (is.call(cs[[2]])) tryCatch(eval(cs[[2]]), error = function(e) e) else cs[[2]]
  r <- list()
  run <- function(expr) tryCatch(withCallingHandlers(expr,
      warning = function(w) invokeRestart("muffleWarning"),
      message = function(m) invokeRestart("muffleMessage")), error = function(e) e)
  if (inherits(fo, "error") || inherits(fam, "error")) {
    out <- if (inherits(fo, "error")) fo else fam
  } else if (identical(ARM, "brms")) {
    out <- run(if (is.null(fam)) brms::make_standata(fo, data = d) else brms::make_standata(fo, data = d, family = fam))
    gp <- run(if (is.null(fam)) brms::get_prior(fo, data = d) else brms::get_prior(fo, data = d, family = fam))
    r$get_prior <- if (inherits(gp, "error")) "error" else "ok"
  } else {
    out <- run(if (is.null(fam)) frm(fo, data = d) else frm(fo, data = d, family = fam))
    gp <- run(if (is.null(fam)) get_prior(fo, data = d) else get_prior(fo, data = d, family = fam))
    r$get_prior <- if (inherits(gp, "error")) paste("error:", substr(conditionMessage(gp), 1, 60)) else "ok"
    if (!inherits(out, "error")) {
      r$ll <- as.numeric(logLik(out)); r$est <- unlist(out$estimates)
      nd <- d[1:5, setdiff(names(d), "n")]
      p1 <- run(predict(out, newdata = nd))
      r$predict_newdata_no_trials <- if (inherits(p1, "error")) paste("error:", substr(conditionMessage(p1), 1, 70)) else "ok"
      p2 <- run(predict(out, newdata = d[1:5, ]))
      r$predict_newdata <- if (inherits(p2, "error")) paste("error:", substr(conditionMessage(p2), 1, 70)) else "ok"
      s1 <- run(simulate(out, nsim = 1, seed = 1))
      r$simulate <- if (inherits(s1, "error")) paste("error:", substr(conditionMessage(s1), 1, 70)) else "ok"
    }
    if (!is.null(fam) && !inherits(fam, "error")) {
      tpl <- run(par_template(fo, data = d, family = fam))
      fs <- if (inherits(tpl, "error")) tpl else run(frm_simulate(fo, data = d, family = fam, newparams = tpl, seed = 1))
      r$frm_simulate <- if (inherits(fs, "error")) paste("error:", substr(conditionMessage(fs), 1, 70)) else "ok"
    }
  }
  r$status <- if (inherits(out, "error")) "error" else "ok"
  r$message <- if (inherits(out, "error")) conditionMessage(out) else ""
  res[[nm]] <- r
  cat(sprintf("%-28s %-5s gp %-6s %s\n", nm, r$status, substr(r$get_prior, 1, 6), substr(r$message, 1, 90)))
}
saveRDS(res, sprintf("dev/famlink-rev2-trials-%s.rds", ARM))
