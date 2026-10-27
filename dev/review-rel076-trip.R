# Review: the record's open item. dev/rel076-overlap-trip.R compares
# rstan::log_prob() of the defaults' stanfit and the table's stanfit at
# the same points, after BOTH runs. tmbstan's sampling method sets
# GLOBAL pointers to the objective (.Call("set_pointers", ...)), so
# after the second run both stanfits evaluate the second run's
# objective, and that comparison cannot see a difference. Here each
# run's objective is evaluated right after its own run, at one fixed
# set of points, and the two are compared bitwise.
#   Rscript dev/review-rel076-trip.R > dev/review-rel076-out/trip.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
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
# points: unconstrained, fixed before either run
set.seed(11)
U <- NULL
eval_at <- function(sf, U) {
  t(apply(U, 1, function(u) {
    c(lp = rstan::log_prob(sf, u), rstan::grad_log_prob(sf, u))
  }))
}
a <- run(NULL)

ra <- rstan::extract(a$stanfit, permuted = FALSE)[, 1, ]
pn <- setdiff(dimnames(ra)[[2]], "lp__")
u0 <- rstan::unconstrain_pars(a$stanfit, list(y = unname(ra[1, pn])))
np <- length(u0)
U <- rbind(t(apply(ra[, pn], 1, function(v) {
  rstan::unconstrain_pars(a$stanfit, list(y = unname(v)))
})), t(vapply(1:300, function(i) u0 + stats::rnorm(np, 0, 1), numeric(np))))
ea <- eval_at(a$stanfit, U)
# a's lp__ at its own draws, recorded during sampling, against the
# evaluation now: the pointer is a's here
cat("a: recorded lp__ equals log_prob now at its 200 draws:",
    isTRUE(all.equal(unname(ra[, "lp__"]), unname(ea[1:200, "lp"]),
                     tolerance = 0)), "\n")
b <- run(dp)
eb <- eval_at(b$stanfit, U)
# and a's stanfit now, after b: if the pointer is global it reads b's
ea_after <- eval_at(a$stanfit, U)
cat("points:", nrow(U), "\n")
cat("a's stanfit after b's run equals b's evaluation (global pointer):",
    identical(ea_after, eb), "\n")
dl <- eb[, "lp"] - ea[, "lp"]
cat("log density, each run's own objective: differs at",
    sum(dl != 0), "of", nrow(U), "| max abs", format(max(abs(dl))),
    "| max relative", format(max(abs(dl) / abs(ea[, "lp"]))), "\n")
dg <- abs(eb[, -1] - ea[, -1])
cat("gradient: differs at", sum(apply(dg, 1, max) != 0), "of", nrow(U),
    "| max abs", format(max(dg)), "\n")
spa <- rstan::get_sampler_params(a$stanfit, inc_warmup = TRUE)[[1]]
spb <- rstan::get_sampler_params(b$stanfit, inc_warmup = TRUE)[[1]]
wa <- rstan::extract(a$stanfit, permuted = FALSE, inc_warmup = TRUE)[, 1, ]
wb <- rstan::extract(b$stanfit, permuted = FALSE, inc_warmup = TRUE)[, 1, ]
dd <- apply(abs(wa - wb), 1, max)
cat("warmup iteration of the first draw difference: above 0 at",
    which(dd > 0)[1], "; above 1e-8 at", which(dd > 1e-8)[1],
    "; above 1e-2 at", which(dd > 1e-2)[1], "of", nrow(wa), "\n")
cat("differences, iterations 1 to 12:", format(dd[1:12], digits = 3), "\n")
cat("posterior means a:", format(colMeans(ra), digits = 5), "\n")
rb <- rstan::extract(b$stanfit, permuted = FALSE)[, 1, ]
cat("posterior means b:", format(colMeans(rb), digits = 5), "\n")
cat("posterior sd a:   ", format(apply(ra, 2, sd), digits = 5), "\n")
# the same objective run twice: reproduces (the record's control)
b2 <- run(dp)
cat("table run twice identical:",
    identical(unname(as.matrix(b2$draws)), unname(as.matrix(b$draws))), "\n")
