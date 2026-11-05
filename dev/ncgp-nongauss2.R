# dev/ncgp-nongauss.R with chains started at frm()'s z as well (x is sorted here, so z is brms's order): Poisson at a large
# log mean, binomial at many trials, bernoulli. frm() (Laplace) against brms's
# posterior, unedited program, brms's default priors on both sides, seed in the
# row. Also the Laplace value against the same fit's z mode by hand. Output
# NCGPZ lines. On the base build gives the frm() row only.
#   Rscript dev/ncgp-nongauss.R <case> <lib or base> [iter]
args <- commandArgs(TRUE)
case <- args[1]
libs <- c("C:/Users/adf44/source/r/rellib-r16",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (args[2] != "base") libs <- c(args[2], libs)
.libPaths(libs)
iter <- if (length(args) > 2) as.integer(args[3]) else 1500L
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
n <- 40L
set.seed(5601)
x <- sort(round(stats::runif(n, 0, 6), 3))
f <- sin(2 * x)
spec <- switch(
  case,
  pois4 = list(y = stats::rpois(n, exp(4 + f)), fam = "poisson"),
  pois12 = list(y = stats::rpois(n, exp(12 + f)), fam = "poisson"),
  pois16 = list(y = stats::rpois(n, exp(16 + f)), fam = "poisson"),
  binom1e7 = list(y = stats::rbinom(n, 1e7, plogis(f)), fam = "binomial", N = 1e7),
  binom1e5 = list(y = stats::rbinom(n, 1e5, plogis(f)), fam = "binomial", N = 1e5),
  bern = list(y = stats::rbinom(200L, 1L, plogis(rep(f, 5L))), fam = "bernoulli"))
if (case == "bern") x <- rep(x, 5L)
d <- data.frame(x = x, y = spec$y)
if (!is.null(spec$N)) d$N <- spec$N
form <- function(ns) {
  bfn <- get("bf", ns)
  if (!is.null(spec$N)) bfn(y | trials(N) ~ gp(x)) else bfn(y ~ gp(x))
}
famf <- function() switch(spec$fam, poisson = poisson(), binomial = binomial(), bernoulli = bernoulli())
t0 <- Sys.time()
w <- character()
fit <- withCallingHandlers(
  frm(form(asNamespace("frmtmb")), data = d, family = famf()),
  warning = function(cnd) { w <<- c(w, substr(conditionMessage(cnd), 1, 70)); invokeRestart("muffleWarning") })
bk <- fit$frame$re_blocks[[1L]]
th <- fit$estimates$theta[bk$theta_idx]
sum0 <- c(b0 = fit$estimates$beta[[1L]], sdgp = exp(th[[1L]]),
          lscale = exp(th[[2L]]) / bk$gp_lscale_div)
cat(sprintf("NCGPZ %s frm ML ll %.4f b0 %.4g sdgp %.4g lscale %.4g conv %d %.1fs warn %s\n",
            case, -fit$opt$objective, sum0[1], sum0[2], sum0[3], fit$opt$convergence,
            as.numeric(Sys.time() - t0, units = "secs"), paste(unique(w), collapse = " | ")))
if (args[2] == "base") quit("no")
suppressPackageStartupMessages(library(brms))
t0 <- Sys.time()
z0 <- as.numeric(fit$estimates$b[bk$b_idx])
b <- suppressWarnings(brms::brm(form(asNamespace("brms")), data = d, family = famf(),
  chains = 4L, cores = 4L, iter = iter, warmup = iter %/% 2L, seed = 5601,
  control = list(adapt_delta = 0.95, max_treedepth = 12L),
  init = function() list(sdgp_1 = array(sum0[["sdgp"]]),
                         lscale_1 = array(sum0[["lscale"]], c(1L, 1L)),
                         Intercept = sum0[["b0"]], zgp_1 = array(z0, length(z0))),
  refresh = 0, silent = 2L))
sm <- as.data.frame(posterior::summarise_draws(posterior::subset_draws(
  posterior::as_draws_df(b), variable = c("sdgp_gpx", "lscale_gpx")),
  ~quantile(.x, c(.05, .5, .95)), "rhat", "ess_bulk"))
print(sm)
np <- brms::nuts_params(b)
cat(sprintf("NCGPZ %s brms %.0fs divergent %d maxdepth %d\n", case,
            as.numeric(Sys.time() - t0, units = "secs"),
            sum(np$Value[np$Parameter == "divergent__"]),
            sum(np$Value[np$Parameter == "treedepth__"] >= 12)))
pr <- frmtmb::default_prior(form(asNamespace("frmtmb")), data = d, family = famf(), route = "sample")
fm <- suppressWarnings(frm(form(asNamespace("frmtmb")), data = d, family = famf(), prior = pr))
th2 <- fm$estimates$theta[bk$theta_idx]
cat(sprintf("NCGPZ %s frm MAP sdgp %.4g lscale %.4g conv %d\n", case, exp(th2[[1L]]),
            exp(th2[[2L]]) / bk$gp_lscale_div, fm$opt$convergence))
lp <- brms::posterior_linpred(b)
pm <- colMeans(lp)
ps <- apply(lp, 2, stats::sd)
eta <- frmtmb::frm_linpred(fm, type = "link")
se <- frmtmb::frm_linpred(fm, type = "link", se.fit = TRUE)$se.fit
z <- (eta - pm) / ps
cat(sprintf(paste0("NCGPZ %s frm MAP eta against posterior linpred: max |z| %.3f, ",
                   "median |z| %.3f; se / posterior sd: median %.3f, range %.3f %.3f\n"),
            case, max(abs(z)), stats::median(abs(z)), stats::median(se / ps),
            min(se / ps), max(se / ps)))
