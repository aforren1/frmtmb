# dev/rel076-overlap.R arms 2 and 3 found the draws under the sampling
# table passed back as prior = not identical to the defaults' draws.
# acprior's test pins the per-entry densities as identical, and says the
# two stacks sum their entries in another order (an ulp apart), which a
# chain amplifies. So: how far apart are the first draws, the
# objective at the same points, and the posterior means.
#   Rscript dev/rel076-overlap-trip.R <lib>
#
# INVALID as a comparison of the two objectives (the release review,
# dev/reviews/2026-10-27-release.md, finding 2): tmbstan holds the
# objective in process-global pointers, so every log_prob() and
# grad_log_prob() below, called after both runs, evaluates the second
# run's objective (dev/rtmb-pitfalls.md item 23). The reviewer's
# dev/review-rel076-trip.R evaluates each right after its own run. Kept
# unchanged otherwise, as the record of what was measured.
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(20261027)
d <- expand.grid(t = 1:6, g = 1:8)
d$x <- stats::rnorm(nrow(d))
e <- unlist(lapply(1:8, function(g) {
  as.numeric(stats::arima.sim(list(ar = 0.9), 6, sd = 0.45))
}))
d$y <- 1 + 0.5 * d$x + e
fo <- y ~ x + ar(t, g, cov = TRUE)
dp <- default_prior(fo, data = d, route = "sample")
run <- function(prior) {
  args <- list(fo, data = d, chains = 1, iter = 400, seed = 3,
               refresh = 0)
  if (!is.null(prior)) args$prior <- prior
  suppressWarnings(suppressMessages(do.call(frm_sample, args)))
}
a <- run(NULL)
b <- run(dp)
ma <- as.matrix(a$draws)
mb <- as.matrix(b$draws)
dif <- apply(abs(ma - mb), 1, max)
cat("draw-by-draw max abs difference, first 10:",
    format(dif[1:10], digits = 3), "\n")
cat("first draw with a difference above 1e-8:",
    which(dif > 1e-8)[1], "of", nrow(ma), "\n")
cat("posterior means, defaults:", format(colMeans(ma), digits = 5), "\n")
cat("posterior means, table:   ", format(colMeans(mb), digits = 5), "\n")
lp <- cbind(defaults = ma[1:5, "lp__"], table = mb[1:5, "lp__"])
print(lp, digits = 15)
# the log density of each run's own Stan model at the same points: the
# first 5 draws of the defaults' run, unconstrained by each model
ra <- rstan::extract(a$stanfit, permuted = FALSE)[, 1, ]
pn <- setdiff(dimnames(ra)[[2]], "lp__")
lpab <- t(vapply(1:5, function(i) {
  v <- ra[i, pn]
  sk <- list(y = numeric(length(pn)))
  pl <- utils::relist(unname(v), utils::as.relistable(sk))
  ua <- rstan::unconstrain_pars(a$stanfit, pl)
  ub <- rstan::unconstrain_pars(b$stanfit, pl)
  c(a = rstan::log_prob(a$stanfit, ua), b = rstan::log_prob(b$stanfit, ub),
    du = max(abs(ua - ub)))
}, numeric(3)))
print(lpab, digits = 15)
cat("log_prob difference (table - defaults):",
    format(lpab[, "b"] - lpab[, "a"], digits = 15), "\n")
ia <- unlist(rstan::get_inits(a$stanfit))
ib <- unlist(rstan::get_inits(b$stanfit))
cat("chain starts: defaults", format(ia, digits = 6), "\n")
cat("chain starts: table   ", format(ib, digits = 6), "\n")
gd <- vapply(1:5, function(i) {
  pl <- list(y = unname(ra[i, pn]))
  ua <- rstan::unconstrain_pars(a$stanfit, pl)
  max(abs(rstan::grad_log_prob(a$stanfit, ua) -
            rstan::grad_log_prob(b$stanfit, ua)))
}, 0)
cat("gradient max abs difference at the 5 points:", format(gd), "\n")
cat("args: defaults", paste(names(a$stanfit@stan_args[[1]]), collapse = ","),
    "\n")
sa <- a$stanfit@stan_args[[1]]; sb <- b$stanfit@stan_args[[1]]
for (nm in union(names(sa), names(sb))) {
  if (!identical(sa[[nm]], sb[[nm]])) cat("stan_args differ:", nm, "\n")
}
cat("lower/upper: defaults", format(a$stanfit@par_dims$y), "\n")
rb <- rstan::extract(b$stanfit, permuted = FALSE)[, 1, ]
pts <- rbind(ra[, pn], rb[, pn])
cmp <- t(apply(pts, 1, function(v) {
  ua <- rstan::unconstrain_pars(a$stanfit, list(y = unname(v)))
  la <- rstan::log_prob(a$stanfit, ua)
  lb <- rstan::log_prob(b$stanfit, ua)
  ga <- rstan::grad_log_prob(a$stanfit, ua)
  gb <- rstan::grad_log_prob(b$stanfit, ua)
  c(dl = lb - la, dg = max(abs(gb - ga)))
}))
cat("points:", nrow(cmp), "| log_prob differs at", sum(cmp[, "dl"] != 0),
    "| max abs", format(max(abs(cmp[, "dl"]))), "| gradient differs at",
    sum(cmp[, "dg"] != 0), "| max abs", format(max(cmp[, "dg"])), "\n")
# the control: the defaults' run twice in this process
a2 <- run(NULL)
cat("defaults run twice, draws identical:",
    identical(unname(as.matrix(a2$draws)), unname(ma)), "\n")
b2 <- run(dp)
cat("table run twice, draws identical:",
    identical(unname(as.matrix(b2$draws)), unname(mb)), "\n")
spa <- rstan::get_sampler_params(a$stanfit, inc_warmup = FALSE)[[1]]
spb <- rstan::get_sampler_params(b$stanfit, inc_warmup = FALSE)[[1]]
cat("stepsize: defaults", format(spa[1, "stepsize__"], digits = 10),
    "table", format(spb[1, "stepsize__"], digits = 10), "\n")
cat("adaptation (defaults):", rstan::get_adaptation_info(a$stanfit)[[1]], "\n")
cat("adaptation (table):   ", rstan::get_adaptation_info(b$stanfit)[[1]], "\n")
cat("dims:", identical(a$stanfit@par_dims, b$stanfit@par_dims), "\n")
cat("obj par: ", format(a$fit$obj$par), " | ", format(b$fit$obj$par), "\n")
cat("identical stanmodel:", identical(a$stanfit@stanmodel@model_code,
                                      b$stanfit@stanmodel@model_code), "\n")
set.seed(11)
u0 <- rstan::unconstrain_pars(a$stanfit, list(y = unname(ra[1, pn])))
wide <- t(vapply(1:400, function(i) {
  u <- u0 + stats::rnorm(length(u0), 0, 3)
  la <- tryCatch(rstan::log_prob(a$stanfit, u), error = function(e) NA)
  lb <- tryCatch(rstan::log_prob(b$stanfit, u), error = function(e) NA)
  c(dl = lb - la, la = la)
}, numeric(2)))
cat("wide points: log_prob differs at", sum(wide[, "dl"] != 0, na.rm = TRUE),
    "of", sum(!is.na(wide[, "dl"])), "| max abs",
    format(max(abs(wide[, "dl"]), na.rm = TRUE)), "\n")
print(head(wide[wide[, "dl"] != 0 & !is.na(wide[, "dl"]), , drop = FALSE]))
