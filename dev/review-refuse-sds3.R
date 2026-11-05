# Item 5: class sds addressed by dpar, nlpar and resp: the taped gap names which
# smooth's sd the prior reached. Also frm_sample's table and prior_draws on the
# same addresses. Seed 9920.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(9920)
n <- 200
d <- data.frame(x = runif(n), z = runif(n), u = runif(n))
d$y <- sin(6 * d$x) + rnorm(n, 0, exp(-1 + d$z))
d$y2 <- cos(6 * d$u) + rnorm(n, 0, 0.5)
gap <- function(lab, form, prior) {
  f0 <- frm(form, data = d)
  fp <- frm(form, data = d, prior = prior)
  par <- f0$opt$par
  g <- fp$obj$fn(par) - f0$obj$fn(par)
  th <- par[names(par) == "theta"]
  labs <- vapply(f0$frame$re_blocks, function(b) paste0(b$dpar %||% "", ":", b$group_name), "")
  each <- vapply(f0$frame$re_blocks, function(b) {
    i <- b$theta_idx[frmtmb:::block_sd_idx(b)]
    -sum(dexp(exp(th[i]), 2, log = TRUE) + th[i])
  }, 0)
  hit <- which(abs(each - g) < 1e-9)
  pr("%-26s gap %.6f; blocks [%s], each %s; reached: %s", lab, g, paste(labs, collapse = " "), paste(round(each, 4), collapse = " "),
     if (length(hit)) paste(labs[hit], collapse = ",") else if (abs(sum(each) - g) < 1e-9) "all" else "NONE MATCH")
}
ex <- function(...) set_prior("exponential(2)", class = "sds", ...)
gap("dpar sigma", bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5)), ex(dpar = "sigma"))
gap("mu by no dpar", bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5)), ex())
gap("coef in sigma", bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5)), ex(dpar = "sigma", coef = "s(z, k = 5)"))
gap("nlpar a", bf(y ~ a + b, a ~ s(x, k = 5), b ~ s(z, k = 5), nl = TRUE), ex(nlpar = "a"))
gap("nlpar coef", bf(y ~ a + b, a ~ s(x, k = 5), b ~ s(z, k = 5), nl = TRUE), ex(nlpar = "b", coef = "s(z, k = 5)"))
gap("resp y2", mvbf(bf(y ~ s(x, k = 5)), bf(y2 ~ s(u, k = 5)), rescor = FALSE), ex(resp = "y2"))
gap("mu, a by-smooth in sigma", bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5)), ex(dpar = "") )

# on the sampling route: which default and which prior reaches what
s <- tr(suppressMessages(suppressWarnings(frm_sample(bf(y ~ s(x, k = 5), sigma ~ s(z, k = 5)), data = d, prior = ex(dpar = "sigma"),
   chains = 1, iter = 300, seed = 2, refresh = 0, sample_prior = "yes"))))
pr("frm_sample dpar sigma: %s", msg(s, 200))
if (!is_err(s)) {
  pd <- prior_draws(s)
  pr("  prior_draws: %s", paste(names(pd), collapse = " "))
  nm <- grep("sds", names(pd), value = TRUE)
  pr("  prior means: %s", paste(sprintf("%s %.3f", nm, colMeans(pd[nm])), collapse = "; "))
  print(as.data.frame(prior_summary(s))[, c("prior", "class", "coef", "dpar")])
}
s <- tr(suppressMessages(suppressWarnings(frm_sample(bf(y ~ a + b, a ~ s(x, k = 5), b ~ 1, nl = TRUE), data = d,
   prior = ex(nlpar = "a") + set_prior("normal(0, 2)", nlpar = "b"), chains = 1, iter = 300, seed = 2, refresh = 0))))
pr("frm_sample nlpar a: %s", msg(s, 200))
if (!is_err(s)) pr("  draws: %s", paste(grep("^sd", names(posterior::as_draws_df(s)), value = TRUE), collapse = " "))
