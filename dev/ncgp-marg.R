# frm() against the posterior of brms's own model on low noise, with converged
# chains. brms's program does not mix at this noise (Rhat 2.78 at 20 rows, 3000
# iterations, tree depth 14), so the field is integrated out, which leaves the same
# posterior of sigma, sdgp and lscale under brms's default priors. The nugget is
# brms's absolute 1e-12 in one arm and frmtmb's relative one in the other. Seed 6601.
# frm() on the base build for the failing reading.
#   Rscript dev/ncgp-marg.R <lib or base>
a <- commandArgs(TRUE)
libs <- c("C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (a[1] != "base") libs <- c(a[1], libs)
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"))
suppressPackageStartupMessages(library(frmtmb))
H <- new.env(parent = asNamespace("frmtmb"))
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
set.seed(6601)
n <- 40
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
mod <- if (a[1] != "base") H$brms_stan_model(code)
for (s in c(3e-3, 1e-3)) {
  d <- data.frame(x = x, y = sin(2 * x) + s * e0)
  gp <- brms::get_prior(brms::bf(y ~ gp(x)), data = d, family = gaussian())
  lp <- gp$prior[gp$class == "lscale" & nzchar(gp$prior)][1]
  pr <- as.numeric(strsplit(gsub("inv_gamma\\(|\\)", "", lp), ",")[[1]])
  ic <- gp$prior[gp$class == "Intercept"][1]
  loc <- as.numeric(strsplit(gsub("student_t\\(|\\)", "", ic), ",")[[1]])[2]
  sd0 <- brms::make_standata(brms::bf(y ~ gp(x)), data = d, family = gaussian())
  xs <- x / sd0$dmax_1
  # brms's priors on frm(), through the sampling route's table, on the base and here
  bp <- c(frmtmb::set_prior("student_t(3, 0, 2.5)", class = "sdgp"),
          frmtmb::set_prior("student_t(3, 0, 2.5)", class = "sigma"),
          frmtmb::set_prior(lp, class = "lscale", coef = "gpx"))
  map <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian(), prior = bp))
  ml <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
  th <- exp(map$estimates$theta)
  cat(sprintf("NCGPMG s %.0e %s frm ML sigma %.4g; MAP sigma %.4g sdgp %.4g lscale %.4g\n", s, a[1], stats::sigma(ml)[1],
              stats::sigma(map)[1], th[1], th[2] / sd0$dmax_1))
  if (is.null(mod)) next
  dat <- list(N = n, y = d$y, x = as.array(xs), ic_loc = loc, shape = pr[1], rate = pr[2])
  for (arm in c("brms", "rel")) {
    dat$rel <- if (arm == "rel") 1e-12 else 0
    dat$abso <- if (arm == "rel") 0 else 1e-12
    fit <- suppressWarnings(rstan::sampling(
      mod, data = dat, chains = 4, cores = 2, iter = 2000, seed = 77, refresh = 0,
      control = list(adapt_delta = 0.95),
      init = function() list(b0 = 0, sdgp = 1, lscale = 0.2, sigma = s)))
    sm <- rstan::summary(fit)$summary
    q <- function(nm) signif(as.numeric(stats::quantile(rstan::extract(fit, nm)[[1]],
                                                        c(0.05, 0.5, 0.95))), 4)
    # the marginal mode, by nlminb on the same density at this arm's nugget
    cat(sprintf("NCGPMG s %.0e %-4s sigma %s (rhat %.3f) sdgp %s (rhat %.3f) lscale %s (rhat %.3f) divergent %d\n",
                s, arm, paste(q("sigma"), collapse = " "), sm["sigma", "Rhat"],
                paste(q("sdgp"), collapse = " "), sm["sdgp", "Rhat"],
                paste(q("lscale"), collapse = " "), sm["lscale", "Rhat"],
                sum(vapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                           function(p) sum(p[, "divergent__"]), 0))))
  }
}
