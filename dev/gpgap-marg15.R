# The target the drawn runs of dev/gpgap-sampler.R aim at: the marginal posterior,
# field integrated out so the chains mix (the review's dev/review-gpgap-marg.R), on
# the same data, seed 6601, at brms's default priors, three nuggets.
#   Rscript dev/gpgap-marg15.R <n> <s>
a <- commandArgs(trailingOnly = TRUE)
n <- as.integer(a[1]); s <- as.numeric(a[2])
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
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
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + s * e0)
gp <- brms::get_prior(brms::bf(y ~ gp(x)), data = d, family = gaussian())
lp <- gp$prior[gp$class == "lscale" & nzchar(gp$prior)][1]
pr <- as.numeric(strsplit(gsub("inv_gamma[(]|[)]", "", lp), ",")[[1]])
ic <- gp$prior[gp$class == "Intercept"][1]
loc <- as.numeric(strsplit(gsub("student_t[(]|[)]", "", ic), ",")[[1]])[2]
sd0 <- brms::make_standata(brms::bf(y ~ gp(x)), data = d, family = gaussian())
xs <- x / sd0$dmax_1
dat <- list(N = n, y = d$y, x = as.array(xs), ic_loc = loc, shape = pr[1], rate = pr[2])
arms <- list("frm_sample before, rel 1e-6" = c(1e-6, 0),
             "frm_sample now, rel 1e-12" = c(1e-12, 0),
             "brms, abs 1e-12" = c(0, 1e-12))
for (nm in names(arms)) {
  dat$rel <- arms[[nm]][1]
  dat$abso <- arms[[nm]][2]
  fit <- suppressWarnings(rstan::sampling(
    mod, data = dat, chains = 4, iter = 4000, seed = 77, refresh = 0, cores = 1,
    control = list(adapt_delta = 0.95),
    init = function() list(b0 = 0, sdgp = 1, lscale = 0.2, sigma = s)))
  sm <- rstan::summary(fit)$summary
  q <- function(v) signif(as.numeric(stats::quantile(rstan::extract(fit, v)[[1]],
                                                     c(0.05, 0.5, 0.95))), 4)
  cat(sprintf("MARG n %d s %.0e %-28s sigma %s (rhat %.3f) sdgp %s (rhat %.3f) divergent %d\n",
              n, s, nm, paste(q("sigma"), collapse = " "), sm["sigma", "Rhat"],
              paste(q("sdgp"), collapse = " "), sm["sdgp", "Rhat"],
              sum(vapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                         function(p) sum(p[, "divergent__"]), 0))))
}
