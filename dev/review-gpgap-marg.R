# Review of lane gpgap, item 5: the posterior of sigma on the sampling route,
# with the field integrated out so the chains mix. frm_sample's density has the
# relative nugget and brms's the absolute one; their marginals differ by that
# one term. Same data (seed 6601) and brms's default priors. The drawn runs,
# which mix badly at this noise, are dev/review-gpgap-sample*.R.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-marg.R marg.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
H <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-stan.R", envir = H)
sys.source("tests/testthat/helper-brms.R", envir = H)
code <- "
data {
  int N; vector[N] y; array[N] real x;
  real rel; real abso;
  real ic_loc; real shape; real rate;
}
parameters {
  real b0; real<lower=0> sdgp; real<lower=0> lscale; real<lower=0> sigma;
}
model {
  matrix[N, N] V = gp_exp_quad_cov(x, sdgp, lscale);
  for (n in 1:N) V[n, n] += rel * square(sdgp) + abso + square(sigma);
  y ~ multi_normal(rep_vector(b0, N), V);
  b0 ~ student_t(3, ic_loc, 2.5);
  sdgp ~ student_t(3, 0, 2.5);
  sigma ~ student_t(3, 0, 2.5);
  lscale ~ inv_gamma(shape, rate);
}
"
mod <- H$brms_stan_model(code)
set.seed(6601)
n <- 40
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
for (s in c(1e-2, 3e-3, 1e-3, 3e-4)) {
  d <- data.frame(x = x, y = sin(2 * x) + s * e0)
  gp <- brms::get_prior(brms::bf(y ~ gp(x)), data = d, family = gaussian())
  lp <- gp$prior[gp$class == "lscale" & nzchar(gp$prior)][1]
  pr <- as.numeric(strsplit(gsub("inv_gamma\\(|\\)", "", lp), ",")[[1]])
  ic <- gp$prior[gp$class == "Intercept"][1]
  loc <- as.numeric(strsplit(gsub("student_t\\(|\\)", "", ic), ",")[[1]])[2]
  sd0 <- brms::make_standata(brms::bf(y ~ gp(x)), data = d, family = gaussian())
  xs <- x / sd0$dmax_1
  dat <- list(N = n, y = d$y, x = as.array(xs), ic_loc = loc, shape = pr[1], rate = pr[2])
  for (arm in c("frm", "brms")) {
    dat$rel <- if (arm == "frm") 1e-6 else 0
    dat$abso <- if (arm == "frm") 0 else 1e-12
    set.seed(1)
    fit <- suppressWarnings(rstan::sampling(
      mod, data = dat, chains = 2, iter = 2000, seed = 77, refresh = 0, cores = 1,
      control = list(adapt_delta = 0.95),
      init = function() list(b0 = 0, sdgp = 1, lscale = 0.2, sigma = s)))
    sm <- rstan::summary(fit)$summary
    q <- function(nm) signif(as.numeric(stats::quantile(rstan::extract(fit, nm)[[1]],
                                                        c(0.05, 0.5, 0.95))), 4)
    cat(sprintf("MARG s %.0e %-4s sigma %s (rhat %.3f) sdgp %s (rhat %.3f) lscale %s divergent %d\n",
                s, arm, paste(q("sigma"), collapse = " "), sm["sigma", "Rhat"],
                paste(q("sdgp"), collapse = " "), sm["sdgp", "Rhat"],
                paste(q("lscale"), collapse = " "),
                sum(vapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                           function(p) sum(p[, "divergent__"]), 0))))
  }
}
