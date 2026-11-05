# Attack 8 and 2: frm_sample() against brms on the same data, default priors on
# both sides, one design per process, with frm() beside them. Also the methods on
# the draws, the default-prior tables, and variables().
# Usage: Rscript review-latres-sample.R <design> <seed>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
design <- args[[1]]
seed <- as.integer(args[[2]])
suppressPackageStartupMessages({library(frmtmb.sample); library(rstan)})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
designs <- list(
  pois = list(fam = "poisson", G = 10, K = 15, b0 = -0.2, b1 = 0.5, s = 0.5, ar = 0.6),
  bern = list(fam = "bernoulli", G = 20, K = 10, b0 = 0.2, b1 = 0.5, s = 1, ar = 0.6),
  bern_long = list(fam = "bernoulli", G = 5, K = 60, b0 = 0.2, b1 = 0.5, s = 1, ar = 0.8),
  binom3 = list(fam = "binomial", G = 15, K = 10, b0 = 0.2, b1 = 0.5, s = 0.8, ar = 0.6, tr = 3),
  pois_ar9 = list(fam = "poisson", G = 10, K = 15, b0 = -0.2, b1 = 0.5, s = 0.3, ar = 0.9),
  bern_F = list(fam = "bernoulli", G = 20, K = 10, b0 = 0.2, b1 = 0.5, s = 1, ar = 0.6, form = "ar(t, g)"),
  pois_F = list(fam = "poisson", G = 10, K = 15, b0 = -0.2, b1 = 0.5, s = 0.5, ar = 0.6, form = "ar(t, g)"))
D <- designs[[design]]
set.seed(seed)
d <- expand.grid(t = seq_len(D$K), g = seq_len(D$G))
d <- d[order(d$g, d$t), ]
d$x <- rnorm(nrow(d))
e <- numeric(nrow(d))
for (g in seq_len(D$G)) {
  idx <- which(d$g == g)
  e[idx[1]] <- rnorm(1, 0, D$s / sqrt(1 - D$ar^2))
  for (i in idx[-1]) e[i] <- D$ar * e[i - 1] + rnorm(1, 0, D$s)
}
eta <- D$b0 + D$b1 * d$x + e
d$y <- switch(D$fam, poisson = rpois(nrow(d), exp(eta)), bernoulli = rbinom(nrow(d), 1, plogis(eta)),
              binomial = rbinom(nrow(d), D$tr, plogis(eta)))
d$tr <- D$tr %||% 1L
d$g <- factor(d$g)
term <- D$form %||% "ar(t, g, cov = TRUE)"
lhs <- if (D$fam == "binomial") "y | trials(tr)" else "y"
form <- as.formula(paste(lhs, "~ x +", term))
famc <- function() get(D$fam, asNamespace("frmtmb"))()

cat("== frm()\n")
f <- frm(frmtmb::bf(form), family = famc(), data = d)
print(summary(f))

fb <- brms::bf(form)
fam_b <- brms::brmsfamily(D$fam)
cat("== default_prior: brms\n")
pb <- as.data.frame(brms::default_prior(fb, data = d, family = fam_b))[, c("prior", "class", "coef", "lb", "ub")]
print(pb)
cat("== default_prior: frm route = sample\n")
ps <- as.data.frame(default_prior(frmtmb::bf(form), family = famc(), data = d, route = "sample"))
print(ps[, intersect(c("prior", "class", "coef", "lb", "ub"), names(ps))])
cat("== default_prior: frm route = fit\n")
pf <- as.data.frame(default_prior(frmtmb::bf(form), family = famc(), data = d))
print(pf[, intersect(c("prior", "class", "coef", "lb", "ub"), names(pf))])

t0 <- proc.time()[[3]]
dr <- frm_sample(frmtmb::bf(form), family = famc(), data = d, chains = 4, iter = 2000, warmup = 1000,
                 seed = 11, refresh = 0)
cat("frm_sample seconds", proc.time()[[3]] - t0, "\n")
print(prior_summary(dr))
v <- variables(dr)
cat("n variables", length(v), "\n")
print(head(v, 8)); print(tail(v, 8))
sm <- as.data.frame(posterior::summarise_draws(posterior::as_draws_df(dr)))
keep <- c("b_Intercept", "b_x", "ar[1]", "sderr")
print(sm[sm$variable %in% keep, c("variable", "mean", "median", "sd", "rhat", "ess_bulk")], digits = 3)

code <- brms::stancode(fb, data = d, family = fam_b)
sdat <- brms::standata(fb, data = d, family = fam_b)
mod <- hel$brms_stan_model(code)
t0 <- proc.time()[[3]]
sf <- rstan::sampling(mod, data = sdat, chains = 4, iter = 2000, warmup = 1000, seed = 11, refresh = 0)
cat("brms seconds", proc.time()[[3]] - t0, "\n")
print(round(rstan::summary(sf, pars = c("b_Intercept", "b", "ar", "sderr"))$summary[, c("mean", "50%", "sd", "Rhat", "n_eff")], 3))
scaf <- brms::brm(fb, d, family = fam_b, empty = TRUE)
scaf$fit <- sf
bfit <- brms::rename_pars(scaf)
cat("brms variables tail:\n"); print(tail(brms::variables(bfit), 3))
cat("brms err named?\n"); print(head(grep("^err", brms::variables(bfit), value = TRUE), 3))

cat("== methods on draws\n")
ep <- msg_of(posterior_epred(dr))
epb <- brms::posterior_epred(bfit)
if (is.null(ep$err)) {
  cat("epred dim", dim(ep$val), "  max |mean(frm)/mean(brms)-1| over rows", max(abs(colMeans(ep$val) / colMeans(epb) - 1)),
      "  median", median(abs(colMeans(ep$val) / colMeans(epb) - 1)), "\n")
} else cat("epred ERR", ep$err, "\n")
pp <- msg_of(posterior_predict(dr, ndraws = 500))
cat("posterior_predict:", if (!is.null(pp$err)) pp$err else paste(dim(pp$val), collapse = " x "), "\n")
cat("   mean of draws", mean(pp$val), " mean y", mean(d$y), "\n")
ll <- msg_of(log_lik(dr))
llb <- brms::log_lik(bfit)
if (is.null(ll$err)) {
  cat("log_lik dim", dim(ll$val), "  sum of pointwise mean: frm", sum(colMeans(ll$val)), " brms", sum(colMeans(llb)), "\n")
} else cat("log_lik ERR", ll$err, "\n")
lo <- msg_of(loo(dr))
lob <- msg_of(brms::loo(bfit))
if (is.null(lo$err)) {
  cat("loo elpd frm", lo$val$estimates["elpd_loo", 1], "se", lo$val$estimates["elpd_loo", 2], "p_loo", lo$val$estimates["p_loo", 1], " brms", lob$val$estimates["elpd_loo", 1], "\n")
  cat("  warnings:", substr(paste(lo$warn, collapse = " | "), 1, 200), "\n")
} else cat("loo ERR", lo$err, "\n")
if (!grepl("cov = TRUE", term)) {
  cat("newdata refusal: "); print(msg_of(posterior_epred(dr, newdata = d[1:3, ]))$err)
} else {
  nd <- d[d$g %in% levels(d$g)[1:2], ]
  a <- msg_of(posterior_epred(dr, newdata = nd))
  b <- msg_of(brms::posterior_epred(bfit, newdata = nd))
  cat("epred newdata frm", mean(a$val), " brms", mean(b$val), " frm ML fitted(newdata)", mean(fitted(f, newdata = nd)[, 1]), " ML predict(newdata)", mean(predict(f, newdata = nd, ndraws = 2000)[, 1]), "\n")
  cat("posterior_predict newdata mean frm", mean(posterior_predict(dr, newdata = nd, ndraws = 2000)), " brms", mean(brms::posterior_predict(bfit, newdata = nd, ndraws = 2000)), "\n")
}
cat("conditional_effects:", class(msg_of(conditional_effects(dr, "x"))$val), "\n")
ty <- msg_of(frm_sample(frmtmb::bf(form), family = famc(), data = d, prior = set_prior("normal(0, 0.5)", class = "sderr"), chains = 1, iter = 300, warmup = 200, seed = 11, refresh = 0))
cat("prior on sderr accepted: ", if (is.null(ty$err)) "yes" else ty$err, "\n")
cat("MATCH-TABLE\n")
cat("sderr brms median", median(rstan::extract(sf, "sderr")[[1]]), " frm_sample median", sm$median[sm$variable == "sderr"], " frm ML", exp(f$estimates$theta[1]), "\n")
cat("done\n")
