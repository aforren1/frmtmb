# coef() on a fit with no group-level term: brms against frm(), on the same model.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages({library(rstan); library(brms)})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:12))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$y <- rpois(n, exp(.5 + .4 * d$x))
pri <- c(set_prior("", class = "sderr"), set_prior("", class = "Intercept"))
form <- y ~ x + ar(t, g, cov = TRUE)
code <- brms::stancode(form, d, family = brms::brmsfamily("poisson"), prior = pri)
sdat <- brms::standata(form, d, family = brms::brmsfamily("poisson"), prior = pri)
m <- hel$brms_stan_model(code)
sf <- rstan::sampling(m, data = sdat, algorithm = "Fixed_param", chains = 1, iter = 20, warmup = 0, refresh = 0,
                      init = list(list(b = array(.4, 1), Intercept = .5, ar = array(.5, 1), zerr = rep(0, n), sderr = .5)))
scaf <- brm(form, d, family = brms::brmsfamily("poisson"), prior = pri, empty = TRUE)
scaf$fit <- sf
bf1 <- brms::rename_pars(scaf)
cat("brms coef: "); r <- msg_of(stats::coef(bf1)); print(r$err %||% names(r$val))
ce <- msg_of(brms::conditional_effects(bf1, "x"))
cat("brms conditional_effects at x grid head: "); print(round(head(ce$val$x$estimate__, 3), 4))
cat("exp(.5 + .4 x) at the grid: "); print(round(head(exp(.5 + .4 * ce$val$x$x), 3), 4))
