# Area 8: fits and readers on a broad set, every number written with %a, one
# file per arm. FRM_ARM=before or after; compared by dev/review-nlfix-a8cmp.R.
source("dev/review-nlfix-common.R")
out <- sprintf("dev/review-nlfix-a8-%s.tsv", arm)
set.seed(8801)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:12, each = 20)),
                f2 = factor(rep(c("a", "b", "c"), length.out = n)), t = rep(1:20, 12))
d$yn <- 1 + 0.5 * d$x + rnorm(n)
d$yg <- exp(0.2 + 0.3 * d$x + rnorm(n, 0, 0.4))
d$yc <- rpois(n, exp(0.3 + 0.4 * d$x))
d$yb <- rbinom(n, 1, plogis(0.3 * d$x))
d$yt <- rbinom(n, 10, plogis(0.3 * d$x)); d$tr <- 10L
d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
d$yo2 <- cut(0.3 * d$z + rlogis(n), c(-Inf, -0.5, 0.5, Inf), labels = 1:3, ordered_result = TRUE)
d$yh <- ordered(ifelse(runif(n) < 0.2, 0L, as.integer(d$yo)), levels = 0:4)
d$ycat <- factor(sample(c("a", "b", "c"), n, TRUE))
m <- t(rmultinom(n, 10, c(1, 1, 1))); colnames(m) <- c("a", "b", "c"); d$cnt <- m
s <- (m + 0.5) / 11.5; d$sim <- s
d$ybeta <- plogis(0.3 * d$x + rnorm(n))
d$gg <- factor(rep(c("p", "q"), length.out = n))
cases <- list(
  gauss = list(y = bf(yn ~ x), fam = gaussian()),
  gauss_sig = list(y = bf(yn ~ x, sigma ~ z), fam = gaussian()),
  gauss_re = list(y = bf(yn ~ x + (1 | g)), fam = gaussian()),
  student = list(y = bf(yn ~ x), fam = student()),
  gamma_log = list(y = bf(yg ~ x), fam = Gamma("log")),
  gamma_inv = list(y = bf(yg ~ x), fam = Gamma()),
  ig = list(y = bf(yg ~ x), fam = inverse.gaussian()),
  poisson = list(y = bf(yc ~ x), fam = poisson()),
  negbin = list(y = bf(yc ~ x, shape ~ z), fam = negbinomial()),
  zip = list(y = bf(yc ~ x, zi ~ z), fam = zero_inflated_poisson()),
  hurdle = list(y = bf(yc ~ x, hu ~ z), fam = hurdle_poisson()),
  binom = list(y = bf(yt | trials(tr) ~ x), fam = binomial()),
  bern = list(y = bf(yb ~ x), fam = bernoulli()),
  beta = list(y = bf(ybeta ~ x, phi ~ z), fam = Beta()),
  cum0 = list(y = bf(yo ~ 1), fam = cumulative()),
  cum = list(y = bf(yo ~ x), fam = cumulative()),
  cum_noint = list(y = bf(yo ~ 0 + x), fam = cumulative()),
  cum_factor = list(y = bf(yo ~ f2), fam = cumulative()),
  cum_cs = list(y = bf(yo ~ cs(x)), fam = cumulative()),
  cum_re = list(y = bf(yo ~ x + (1 | g)), fam = cumulative()),
  cum_disc = list(y = bf(yo ~ x, disc ~ 0 + z), fam = cumulative()),
  cum_disc_int = list(y = bf(yo ~ x, disc ~ z), fam = cumulative()),
  cum_probit = list(y = bf(yo ~ x), fam = cumulative("probit")),
  cum_equi = list(y = bf(yo ~ x), fam = cumulative(threshold = "equidistant")),
  cum_stz = list(y = bf(yo ~ x), fam = cumulative(threshold = "sum_to_zero")),
  cum_thres = list(y = bf(yo | thres(3) ~ x), fam = cumulative()),
  cum_thres_gr = list(y = bf(yo | thres(gr = gg) ~ x), fam = cumulative()),
  sratio = list(y = bf(yo ~ x), fam = sratio()),
  cratio = list(y = bf(yo ~ x), fam = cratio()),
  acat = list(y = bf(yo ~ x), fam = acat()),
  acat_cs = list(y = bf(yo ~ cs(x)), fam = acat()),
  hurdle_cum = list(y = bf(yh ~ x, hu ~ z), fam = hurdle_cumulative()),
  mix_cum = list(y = bf(yo ~ x), fam = mixture(cumulative(), cumulative())),
  mix_gauss = list(y = bf(yn ~ x), fam = mixture(gaussian(), gaussian())),
  cat = list(y = bf(ycat ~ x), fam = categorical()),
  mn = list(y = bf(cnt | trials(tr) ~ x), fam = multinomial()),
  dir = list(y = bf(sim ~ x), fam = dirichlet()),
  mv_ord = list(y = mvbf(bf(yo ~ x, family = cumulative()), bf(yo2 ~ z, family = cumulative()), rescor = FALSE), fam = NULL),
  mv_ord_int = list(y = mvbf(bf(yo ~ 1, family = cumulative()), bf(yn ~ 1, family = gaussian()), rescor = FALSE), fam = NULL),
  mv_gauss = list(y = mvbf(bf(yn ~ x), bf(yg ~ z), rescor = TRUE), fam = gaussian()),
  ar_cov = list(y = bf(yn ~ x + ar(time = t, gr = g, cov = TRUE)), fam = gaussian()),
  ar_cond = list(y = bf(yn ~ x + ar(time = t, gr = g)), fam = gaussian()),
  cosy = list(y = bf(yn ~ x + cosy(time = t, gr = g)), fam = gaussian()),
  nl_gauss = list(y = bf(yn ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), fam = gaussian()),
  nl_exp = list(y = bf(yn ~ a * exp(b * x), a + b ~ 1, nl = TRUE), fam = gaussian(), start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.1))),
  nl_re = list(y = bf(yn ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE), fam = gaussian()),
  nl_pois = list(y = bf(yc ~ a + b * x, a + b ~ 1, nl = TRUE), fam = poisson()),
  nl_sig = list(y = bf(yn ~ a + b * x, sigma ~ z, a + b ~ 1, nl = TRUE), fam = gaussian()),
  nl_gamma_log = list(y = bf(yg ~ a + b * x, a + b ~ 1, nl = TRUE), fam = Gamma("log")),
  nl_prior = list(y = bf(yn ~ a + b * x, a + b ~ 1, nl = TRUE), fam = gaussian(), prior = set_prior("normal(1, 2)", nlpar = "a")),
  nl_chain = list(y = bf(yn ~ eta * 2) + nlf(eta ~ a + b * x) + lf(a ~ 1) + lf(b ~ 1), fam = gaussian()),
  nlf_sigma = list(y = bf(yn ~ x) + nlf(sigma ~ exp(s0 + s1 * z)) + lf(s0 ~ 1) + lf(s1 ~ 1), fam = gaussian()),
  nl_bernoulli = list(y = bf(yb ~ a + b * x, a + b ~ 1, nl = TRUE), fam = bernoulli())
)
hex <- function(v) paste(sprintf("%a", as.numeric(v)), collapse = ",")
con <- file(out, "w")
wr <- function(id, what, v) writeLines(paste(id, what, hex(v), sep = "\t"), con)
for (id in names(cases)) {
  cs <- cases[[id]]
  args <- c(list(cs$y, data = d), if (!is.null(cs$fam)) list(family = cs$fam),
            if (!is.null(cs$start)) list(start = cs$start), if (!is.null(cs$prior)) list(prior = cs$prior))
  r <- tt(do.call(frm, args))
  if (iserr(r)) {
    writeLines(paste(id, "ERROR", gsub("[\r\n\t]", " ", r$msg), sep = "\t"), con)
    next
  }
  wr(id, "logLik", as.numeric(logLik(r)))
  writeLines(paste(id, "names", paste(unlist(lapply(r$estimates, names)), collapse = ","), sep = "\t"), con)
  wr(id, "est", unlist(r$estimates))
  writeLines(paste(id, "warn", gsub("[\r\n\t]", " ", paste(c(attr(r, "w"), attr(r, "m")), collapse = " | ")), sep = "\t"), con)
  v <- tt(as.matrix(vcov(r))); if (!iserr(v)) wr(id, "vcov", v)
  for (m in c("fitted", "frm_linpred")) {
    v <- tt(get(m)(r, newdata = d[1:6, ]))
    if (iserr(v)) writeLines(paste(id, m, "ERR", gsub("[\r\n\t]", " ", substr(v$msg, 1, 200)), sep = "\t"), con) else wr(id, m, unlist(v))
  }
  v <- tt({ set.seed(3); predict(r, newdata = d[1:6, ]) }); if (iserr(v)) writeLines(paste(id, "predict", "ERR", gsub("[\r\n\t]", " ", substr(v$msg, 1, 200)), sep = "\t"), con) else wr(id, "predict", unlist(v))
  v <- tt({ set.seed(4); s <- simulate(r, nsim = 2); suppressWarnings(as.numeric(unlist(lapply(s, as.numeric)))) }); if (iserr(v)) writeLines(paste(id, "simulate", "ERR", gsub("[\r\n\t]", " ", substr(v$msg, 1, 200)), sep = "\t"), con) else wr(id, "simulate", v)
  p <- tt(get_prior(r)); if (!iserr(p)) writeLines(paste(id, "prior", paste(do.call(paste, as.data.frame(p)[, 1:5]), collapse = ";"), sep = "\t"), con)
  p <- tt(summary(r)$fixed); if (!iserr(p) && is.matrix(p)) writeLines(paste(id, "summary", paste(rownames(p), collapse = ";"), sep = "\t"), con)
}
close(con)
