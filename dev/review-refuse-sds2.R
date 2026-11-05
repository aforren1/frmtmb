# Item 5, the rest: frm_simulate, brms's tables, the exact prior contribution
# per smooth sd on the tape, and frm_sample against brms on the posterior of
# sds. Seed 9910.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(9910)
n <- 200
d <- data.frame(x = runif(n), z = runif(n), u = runif(n),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(sample(1:8, n, TRUE)))
d$y <- sin(6 * d$x) + d$z^2 + rnorm(8, 0, 0.4)[d$g] + rnorm(n, 0, 0.4)
d$y2 <- cos(6 * d$u) + rnorm(n, 0, 0.5)

# the exact contribution on the tape, at one point
tape_gap <- function(lab, form, prior, want_idx = NULL) {
  f0 <- frm(form, data = d)
  fp <- frm(form, data = d, prior = prior)
  par <- f0$opt$par
  # the two objectives at one outer point: the gap is the prior alone
  nth <- which(names(par) == "theta")
  tn <- f0$frame$re_blocks
  gap <- fp$obj$fn(par) - f0$obj$fn(par)
  ths <- par[nth]
  # which theta are smooth sd's, from the block table
  sm <- integer(); rg <- integer()
  for (bk in tn) {
    ii <- bk$theta_idx[frmtmb:::block_sd_idx(bk)]
    if (identical(bk$covstruct, "smooth")) sm <- c(sm, ii) else rg <- c(rg, ii)
  }
  want <- -sum(dexp(exp(ths[sm]), 2, log = TRUE) + ths[sm])
  pr("%-34s smooth sd's %d, others %d; objective gap %.10f; exponential(2) on the smooth sd's with Jacobian %.10f; difference %.2e",
     lab, length(sm), length(rg), gap, want, gap - want)
}
ex2 <- set_prior("exponential(2)", class = "sds")
tape_gap("s(x)", bf(y ~ s(x, k = 5)), ex2)
tape_gap("t2(x, z) class", bf(y ~ t2(x, z, k = 4)), ex2)
tape_gap("t2(x, z) coef", bf(y ~ t2(x, z, k = 4)), set_prior("exponential(2)", class = "sds", coef = "t2(x, z, k = 4)"))
tape_gap("s(x, by = f) class", bf(y ~ f + s(x, by = f, k = 4)), ex2)
tape_gap("s(x, by = f) coef", bf(y ~ f + s(x, by = f, k = 4)), set_prior("exponential(2)", class = "sds", coef = "s(x, by = f, k = 4)"))
tape_gap("s(x) + (1 | g), sds only", bf(y ~ s(x, k = 5) + (1 | g)), ex2)
tape_gap("two smooths, class", bf(y ~ s(x, k = 5) + s(z, k = 5)), ex2)
f0 <- frm(bf(y ~ s(x, k = 5) + (1 | g)), data = d)
fp <- frm(bf(y ~ s(x, k = 5) + (1 | g)), data = d, prior = set_prior("exponential(2)", class = "sd"))
par <- f0$opt$par
gap <- fp$obj$fn(par) - f0$obj$fn(par)
bks <- f0$frame$re_blocks
th <- par[names(par) == "theta"]
ig <- unlist(lapply(bks, function(b) if (!identical(b$covstruct, "smooth")) b$theta_idx))
pr("class sd beside a smooth: gap %.10f, on g's sd only %.10f, diff %.2e", gap, -sum(dexp(exp(th[ig]), 2, log = TRUE) + th[ig]),
   gap + sum(dexp(exp(th[ig]), 2, log = TRUE) + th[ig]))

# frm_simulate
fs <- tr(frm_simulate(bf(y ~ s(x, k = 5)), newdata = d, prior = set_prior("exponential(2)", class = "sds") + set_prior("normal(0, 1)", class = "b") +
                        set_prior("normal(0, 1)", class = "Intercept") + set_prior("exponential(1)", class = "sigma"), nsim = 2, seed = 3))
pr("frm_simulate, sds prior: %s", msg(fs, 240))
fs0 <- tr(frm_simulate(bf(y ~ (1 | g)), newdata = d, prior = set_prior("exponential(2)", class = "sd") +
                        set_prior("normal(0, 1)", class = "Intercept") + set_prior("exponential(1)", class = "sigma"), nsim = 2, seed = 3))
pr("frm_simulate, sd prior on RE (control): %s", msg(fs0, 240))
fs2 <- tr(frm_simulate(bf(y ~ s(x, k = 5)), newdata = d, prior = set_prior("exponential(2)", class = "sd") + set_prior("normal(0, 1)", class = "b") +
                        set_prior("normal(0, 1)", class = "Intercept") + set_prior("exponential(1)", class = "sigma"), nsim = 2, seed = 3))
pr("frm_simulate, sd prior on a smooth: %s", msg(fs2, 240))

# brms's own table with a density on sds into frm()
for (fm in c("bf(y ~ t2(x, z, k = 4))", "bf(y ~ f + s(x, by = f, k = 4))", "bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5))")) {
  tb <- brms::default_prior(eval(parse(text = fm), envir = asNamespace("brms")), data = d)
  tb$prior[tb$class == "sds" & !nzchar(tb$coef)] <- "exponential(2)"
  a <- tr(frm(eval(parse(text = fm), envir = asNamespace("frmtmb")), data = d, prior = tb))
  pr("brms table into frm(): %-44s %s", fm, msg(a, 180))
  if (!is_err(a)) {
    a0 <- frm(eval(parse(text = fm), envir = asNamespace("frmtmb")), data = d)
    pr("   objective gap at one point: %.6f", a$obj$fn(a0$opt$par) - a0$obj$fn(a0$opt$par))
  }
}
# the sampling route's table against brms
fm <- bf(y ~ s(x, k = 5) + (1 | g))
pa <- default_prior(fm, data = d, route = "sample")
pb <- brms::default_prior(brms::bf(y ~ s(x, k = 5) + (1 | g)), data = d)
cols <- c("prior", "class", "coef", "group", "lb")
pr("sampling-route rows:"); print(as.data.frame(pa)[pa$class %in% c("sd", "sds"), cols])
pr("brms rows:"); print(as.data.frame(pb)[pb$class %in% c("sd", "sds"), cols])

# posterior of sds on the two samplers, one prior, same data
fb <- brms::bf(y ~ s(x, k = 5) + (1 | g))
pri <- brms::set_prior("exponential(2)", class = "sds") + brms::set_prior("exponential(1)", class = "sd") +
  brms::set_prior("normal(0, 2)", class = "Intercept") + brms::set_prior("exponential(1)", class = "sigma") +
  brms::set_prior("normal(0, 2)", class = "b")
code <- brms::make_stancode(fb, data = d, prior = pri, sample_prior = "yes")
sdat <- E$brms_standata(fb, data = d, prior = pri, sample_prior = "yes")
sb <- suppressWarnings(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 4, iter = 2000, refresh = 0, seed = 11, cores = 1))
eb <- rstan::extract(sb)
pr("brms: names %s", paste(names(eb), collapse = " "))
pfrm <- set_prior("exponential(2)", class = "sds") + set_prior("exponential(1)", class = "sd") +
  set_prior("normal(0, 2)", class = "Intercept") + set_prior("exponential(1)", class = "sigma") +
  set_prior("normal(0, 2)", class = "b")
sm <- frm_sample(bf(y ~ s(x, k = 5) + (1 | g)), data = d, prior = pfrm, chains = 4, iter = 2000, seed = 11, refresh = 0, cores = 1, sample_prior = "yes")
dm <- posterior::as_draws_df(sm)
pr("frm_sample names: %s", paste(grep("^(sd|prior)", names(dm), value = TRUE), collapse = " "))
pd <- prior_draws(sm)
pr("prior_draws names: %s; mean sds %.3f (want 0.5), sd_g %.3f (want 1), sigma %.3f (want 1)", paste(names(pd), collapse = " "),
   mean(pd[[grep("sds", names(pd))[1]]]), mean(pd[[grep("sd_g", names(pd))[1]]]), mean(pd[["sigma"]]))
sb1 <- eb$sds_1
pr("sds posterior mean (sd): brms %.4f (%.4f), frm_sample %.4f (%.4f)", mean(sb1), sd(sb1), mean(dm[["sds_sx_1"]]), sd(dm[["sds_sx_1"]]))
pr("sigma: brms %.4f, frm_sample %.4f", mean(eb$sigma), mean(dm[["sigma"]]))
pr("sd_g: brms %.4f, frm_sample %.4f", mean(eb$sd_1), mean(dm[["sd_g__Intercept"]]))
h <- tr(hypothesis(sm, "sds_sx_1 = 0", class = NULL))
pr("point hypothesis: %s %s", msg(h, 140), if (is_err(h)) "" else format(h$hypothesis$Evid.Ratio, digits = 3))
