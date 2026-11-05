# frm() against brms's posterior on a low-noise design, with chains that
# converge: brms's own program, unedited, brms's default priors on both sides
# (frm() reads them through prior = as a MAP), started at frm()'s field, seed
# 20261105. On 0.78.0 frm()'s sigma reads 25 percent low at true 3e-3 and 3
# orders low at 1e-3. Output NCGPP lines.
#   Rscript dev/ncgp-post.R <true sigma> <n> <lib or base> [iter] [delta]
args <- commandArgs(TRUE)
s <- as.numeric(args[1])
n <- as.integer(args[2])
libs <- c("C:/Users/adf44/source/r/rellib-r16",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (args[3] != "base") libs <- c(args[3], libs)
.libPaths(libs)
iter <- if (length(args) > 3) as.integer(args[4]) else 1500L
delta <- if (length(args) > 4) as.numeric(args[5]) else 0.99
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261105)
x <- round(stats::runif(n, 0, 6), 3)
d <- data.frame(x = x, y = sin(2 * x) + s * stats::rnorm(n))
fit <- suppressWarnings(frm(frmtmb::bf(y ~ gp(x)), data = d, family = gaussian()))
bk <- fit$frame$re_blocks[[1L]]
th <- fit$estimates$theta[bk$theta_idx]
sum0 <- c(sigma = stats::sigma(fit)[1L], sdgp = exp(th[[1L]]), rho = exp(th[[2L]]) / bk$gp_lscale_div)
cat(sprintf("NCGPP s=%.0e n=%d frm ML sigma %.4g sdgp %.4g lscale %.4g conv %d\n",
            s, n, sum0[1], sum0[2], sum0[3], fit$opt$convergence))
if (args[3] == "base") quit("no")
suppressPackageStartupMessages(library(brms))
t0 <- Sys.time()
b <- suppressWarnings(brms::brm(brms::bf(y ~ gp(x)), data = d, family = gaussian(),
  chains = 4L, cores = 4L, iter = iter, warmup = iter %/% 2L, seed = 20261105,
  control = list(adapt_delta = delta, max_treedepth = 14L),
  init = function() list(sdgp_1 = array(sum0[["sdgp"]]),
                         lscale_1 = array(sum0[["rho"]], c(1L, 1L)),
                         sigma = sum0[["sigma"]],
                         Intercept = mean(d$y)),
  refresh = 0, silent = 2L))
tb <- as.numeric(Sys.time() - t0, units = "secs")
sm <- posterior::summarise_draws(posterior::subset_draws(
  posterior::as_draws_df(b), variable = c("sigma", "sdgp_gpx", "lscale_gpx")),
  ~quantile(.x, c(.05, .5, .95)), "rhat", "ess_bulk")
print(as.data.frame(sm))
np <- brms::nuts_params(b)
cat(sprintf("NCGPP s=%.0e n=%d brms %.0fs divergent %d maxdepth %d\n", s, n, tb,
            sum(np$Value[np$Parameter == "divergent__"]),
            sum(np$Value[np$Parameter == "treedepth__"] >= 14)))
# frm()'s MAP at brms's priors: the mode of the marginal of the same model
pr <- frmtmb::default_prior(frmtmb::bf(y ~ gp(x)), data = d, family = gaussian(), route = "sample")
fm <- suppressWarnings(frm(frmtmb::bf(y ~ gp(x)), data = d, family = gaussian(),
                           prior = pr))
th2 <- fm$estimates$theta[bk$theta_idx]
cat(sprintf("NCGPP s=%.0e n=%d frm MAP sigma %.4g sdgp %.4g lscale %.4g conv %d\n",
            s, n, stats::sigma(fm)[1L], exp(th2[[1L]]), exp(th2[[2L]]) / bk$gp_lscale_div,
            fm$opt$convergence))
ep <- brms::posterior_epred(b)
pm <- colMeans(ep)
ps <- apply(ep, 2, stats::sd)
z <- (fitted(fit)[, "Estimate"] - pm) / ps
cat(sprintf("NCGPP s=%.0e n=%d fitted against epred: max |z| %.3f, median |z| %.3f\n",
            s, n, max(abs(z)), stats::median(abs(z))))
ndn <- data.frame(x = c(0.77, 2.11, 4.44))
ep <- brms::posterior_epred(b, newdata = ndn)
z2 <- (frmtmb::frm_linpred(fit, newdata = ndn) - colMeans(ep)) / apply(ep, 2, stats::sd)
cat(sprintf("NCGPP s=%.0e n=%d kriging against epred at new x: |z| %s\n", s, n,
            paste(format(abs(z2), digits = 3), collapse = " ")))
