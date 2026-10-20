# Reviewer of lane sampleprior: the reviewer's own cases (data seed
# 20261114), each a frmtmb call and the brms call it translates.
set.seed(20261114)
n <- 50
x <- rnorm(n, 3, 1)        # mean 3: centering moves the intercept
g <- factor(rep(1:5, length.out = n))
gg <- factor(rep(c("a", "b"), length.out = n))
p01 <- plogis(-0.5 + 0.3 * (x - 3))
D <- list(
  wide = data.frame(x, g, y = 20 + 5 * x + rnorm(n, 0, 10)),
  wpos = data.frame(x, g, y = exp(rnorm(n, 1 + 0.2 * x, 4))),
  lwide = data.frame(x, g, y = exp(rnorm(n, 0.3 * x, 5))),
  cnt = data.frame(x, g, y = rnbinom(n, mu = exp(0.5 + 0.2 * x), size = 2)),
  unit = data.frame(x, g, y = rbeta(n, 20 * p01, 20 * (1 - p01))),
  real = data.frame(x, g, y = 1 + x + rt(n, 4)),
  zcnt = data.frame(x, g, y = ifelse(runif(n) < 0.3, 0,
                                     rpois(n, exp(0.2 * x)))),
  zpos = data.frame(x, g, y = ifelse(runif(n) < 0.25, 0,
                                     exp(rnorm(n, 0.2 * x)))),
  gev = data.frame(x, g, y = 2 + 0.5 * x - log(-log(runif(n)))),
  ord = data.frame(x, g, gg, y = cut(0.7 * x + rlogis(n),
                                     quantile(0.7 * x + rlogis(n),
                                              c(0, .25, .5, .75, 1)) +
                                       c(-Inf, 0, 0, 0, Inf),
                                     labels = FALSE)),
  mix = data.frame(x, g, y = ifelse(runif(n) < 0.4, rnorm(n, 0, 1),
                                    rnorm(n, 5, 1)) + 0.2 * x),
  mv = data.frame(x, g, y1 = 1 + x + rt(n, 5), y2 = 2 - x + rt(n, 5)),
  sqrtc = data.frame(x, g, y = rpois(n, (1 + 0.5 * x)^2)),
  shift = data.frame(x, g, y = 0.3 + exp(rnorm(n, -1 + 0.2 * x, 0.5)))
)
# wiener data from frmtmb.eam's own simulator if present, otherwise a
# crude stand-in with positive times and two choices
D$wien <- data.frame(x, g, rt = 0.35 + rgamma(n, 2, 4),
                     resp = rbinom(n, 1, 0.6))

# frm: list(formula, data, family); brm: list(formula, family)
CASES <- list(
  gauss_wide = list(frm = quote(list(y ~ x, D$wide, gaussian())),
                    brm = quote(list(y ~ x, gaussian()))),
  sigma_x = list(frm = quote(list(bf(y ~ x, sigma ~ x), D$wpos, gaussian())),
                 brm = quote(list(brms::bf(y ~ x, sigma ~ x), gaussian()))),
  lognormal_wide = list(frm = quote(list(y ~ x, D$lwide,
                                         frmtmb::brmsfamily("lognormal"))),
                        brm = quote(list(y ~ x, brms::lognormal()))),
  negbin = list(frm = quote(list(y ~ x, D$cnt,
                                 frmtmb::brmsfamily("negbinomial"))),
                brm = quote(list(y ~ x, brms::negbinomial()))),
  beta_phi_x = list(frm = quote(list(bf(y ~ x, phi ~ x), D$unit,
                                     frmtmb::brmsfamily("beta"))),
                    brm = quote(list(brms::bf(y ~ x, phi ~ x),
                                     brms::Beta()))),
  student = list(frm = quote(list(y ~ x, D$real,
                                  frmtmb::brmsfamily("student"))),
                 brm = quote(list(y ~ x, brms::student()))),
  zip_zi_x = list(frm = quote(list(bf(y ~ x, zi ~ x), D$zcnt,
                         frmtmb::brmsfamily("zero_inflated_poisson"))),
                  brm = quote(list(brms::bf(y ~ x, zi ~ x),
                                   brms::zero_inflated_poisson()))),
  hurdle_ln = list(frm = quote(list(y ~ x, D$zpos,
                                    frmtmb::brmsfamily("hurdle_lognormal"))),
                   brm = quote(list(y ~ x, brms::hurdle_lognormal()))),
  hurdle_ln_hu1 = list(frm = quote(list(bf(y ~ x, hu ~ 1), D$zpos,
                                    frmtmb::brmsfamily("hurdle_lognormal"))),
                   brm = quote(list(brms::bf(y ~ x, hu ~ 1),
                                    brms::hurdle_lognormal()))),
  gev_xi1 = list(frm = quote(list(bf(y ~ x, xi ~ 1), D$gev,
                                  frmtmb::brmsfamily("gen_extreme_value"))),
                 brm = quote(list(brms::bf(y ~ x, xi ~ 1),
                                  brms::gen_extreme_value()))),
  cum_thres_gr = list(frm = quote(list(y | thres(gr = gg) ~ x, D$ord,
                                       frmtmb::brmsfamily("cumulative"))),
                      brm = quote(list(y | thres(gr = gg) ~ x,
                                       brms::cumulative()))),
  mixture = list(frm = quote(list(y ~ x, D$mix,
                                  frmtmb::mixture(gaussian(), gaussian()))),
                 brm = quote(list(y ~ x,
                                  brms::mixture(gaussian(), gaussian())))),
  mv_rescor = list(frm = quote(list(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE),
                                    D$mv, frmtmb::brmsfamily("student"))),
                   brm = quote(list(brms::mvbf(y1 ~ x, y2 ~ x) +
                                      brms::set_rescor(TRUE),
                                    brms::student()))),
  poisson_sqrt = list(frm = quote(list(y ~ x, D$sqrtc,
                                       frmtmb::brmsfamily("poisson",
                                                          link = "sqrt"))),
                      brm = quote(list(y ~ x, brms::brmsfamily("poisson",
                                                         link = "sqrt")))),
  beta_phi_id = list(frm = quote(list(y ~ x, D$unit,
                       frmtmb::brmsfamily("beta", link_phi = "identity"))),
                     brm = quote(list(y ~ x, brms::brmsfamily("beta",
                                               link_phi = "identity")))),
  sigma_id_x = list(frm = quote(list(bf(y ~ x, sigma ~ x), D$wide,
                     frmtmb::brmsfamily("gaussian", link_sigma = "identity"))),
                    brm = quote(list(brms::bf(y ~ x, sigma ~ x),
                      brms::brmsfamily("gaussian",
                                       link_sigma = "identity")))),
  re_sigma = list(frm = quote(list(bf(y ~ x + (1 + x | g), sigma ~ (1 | g)),
                                   D$wide, gaussian())),
                  brm = quote(list(brms::bf(y ~ x + (1 + x | g),
                                            sigma ~ (1 | g)), gaussian()))),
  shifted_ln = list(frm = quote(list(y ~ x, D$shift,
                                     frmtmb::brmsfamily("shifted_lognormal"))),
                    brm = quote(list(y ~ x, brms::shifted_lognormal())))
)

suppressPackageStartupMessages(library(frmtmb.eam))
set.seed(20261115)
D$wien <- frmtmb.eam::ddm_simulate(n, mu = 0.8, bs = 1.4, ndt = 0.3,
                                   bias = 0.5)
D$wien$x <- x
CASES$wiener <- list(frm = quote(list(rt | dec(upper) ~ x, D$wien,
                                      frmtmb.eam::wiener())),
                     brm = quote(list(rt | dec(upper) ~ x, brms::wiener())))
CASES$cum <- list(frm = quote(list(y ~ x, D$ord,
                                   frmtmb::brmsfamily("cumulative"))),
                  brm = quote(list(y ~ x, brms::cumulative())))
CASES$mv_rescor$frm <- quote(list(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE),
                                  D$mv, frmtmb::brmsfamily("student")))
