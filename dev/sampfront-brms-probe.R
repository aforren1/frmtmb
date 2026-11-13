# Lane sampfront: what brms 2.23.0 answers, for the parity items.
#   Rscript dev/sampfront-brms-probe.R > dev/sampfront-out/brms-probe.txt
# Data seed 11; brms seed 1, 1 chain, 400 iterations.
#  A. log_lik(newdata = subset of the fitted rows) against the fitted
#     columns: a plain model, a model with ar(), and loo(newdata =).
#  B. log_lik(combine = FALSE) without rescor: dimensions and dimnames.
#  C. prior_samples(): its warning and its result.
#  D. posterior_average() refusal wording; pp_average(), model_weights().
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
q <- function(...) suppressWarnings(suppressMessages(
  brm(..., chains = 1, iter = 400, refresh = 0, seed = 1,
      backend = "rstan", silent = 2)))
mx <- function(a, b) max(abs(a - b))

set.seed(11)
d <- data.frame(g = factor(rep(1:6, each = 8)), x = rnorm(48),
                t = rep(1:8, 6))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$g] * 0.7 + rnorm(48)
d$z <- -1 + 0.3 * d$x + rnorm(48)

cat("\n== A1 plain y ~ x + (1 | g)\n")
f1 <- q(y ~ x + (1 | g), data = d, prior = prior(normal(0, 5), class = b))
ll <- log_lik(f1)
ls1 <- log_lik(f1, newdata = d[1:10, ])
cat("dim fitted", dim(ll), " dim subset", dim(ls1), "\n")
cat("max |subset - fitted[, 1:10]|", mx(ls1, ll[, 1:10]), "\n")
ls2 <- log_lik(f1, newdata = d[c(5, 3, 3, 20), ])
cat("rows 5,3,3,20: max diff", mx(ls2, ll[, c(5, 3, 3, 20)]), "\n")
l0 <- loo(f1)
l1 <- suppressWarnings(loo(f1, newdata = d[1:10, ]))
cat("loo subset elpd", l1$estimates[1, 1], " sum fitted pointwise 1:10",
    sum(l0$pointwise[1:10, "elpd_loo"]), "\n")
cat("max |pointwise diff|",
    mx(l1$pointwise[, "elpd_loo"], l0$pointwise[1:10, "elpd_loo"]), "\n")
w1 <- suppressWarnings(waic(f1, newdata = d[1:10, ]))
w0 <- suppressWarnings(waic(f1))
cat("waic subset max pointwise diff",
    mx(w1$pointwise[, "elpd_waic"], w0$pointwise[1:10, "elpd_waic"]), "\n")
dw <- d[1:10, ]
dw$x[2] <- dw$x[2] + 1
lw <- log_lik(f1, newdata = dw)
cat("moved x in row 2: max diff row 2", mx(lw[, 2], ll[, 2]),
    " other rows", mx(lw[, -2], ll[, c(1, 3:10)]), "\n")

cat("\n== A2 ar(time = t, gr = g)\n")
f2 <- q(y ~ x + ar(time = t, gr = g), data = d)
ll2 <- log_lik(f2)
sub <- d[d$t >= 3, ]
lsub <- log_lik(f2, newdata = sub)
cat("rows t >= 3: max |subset - fitted cols|",
    mx(lsub, ll2[, which(d$t >= 3)]), "\n")
sub2 <- d[1:10, ]
lsub2 <- log_lik(f2, newdata = sub2)
cat("rows 1:10: max |subset - fitted cols|", mx(lsub2, ll2[, 1:10]), "\n")

cat("\n== B mvbf without rescor, combine = FALSE\n")
f3 <- q(bf(y ~ x) + bf(z ~ x) + set_rescor(FALSE), data = d)
lc <- log_lik(f3, combine = FALSE)
cat("class", class(lc), " dim", dim(lc), "\n")
cat("dimnames:\n"); str(dimnames(lc))
lt <- log_lik(f3)
cat("sum over 3rd margin equals combined:",
    mx(apply(lc, c(1, 2), sum), lt), "\n")
cat("layer y equals resp = 'y':", mx(lc[, , 1], log_lik(f3, resp = "y")),
    "\n")
l1r <- log_lik(f3, resp = "y", combine = FALSE)
cat("one resp combine = FALSE: class", class(l1r), " dim", dim(l1r), "\n")
l2r <- log_lik(f3, resp = c("y", "z"), combine = FALSE)
cat("two resp combine = FALSE: dim", dim(l2r), "\n")
r <- tryCatch(loo(f3, combine = FALSE), error = function(e) e)
cat("loo(combine = FALSE):", if (inherits(r, "error"))
  conditionMessage(r) else class(r), "\n")

cat("\n== C prior_samples()\n")
f4 <- q(y ~ x, data = d, prior = prior(normal(0, 2), class = b),
        sample_prior = "yes")
w <- NULL
ps <- withCallingHandlers(prior_samples(f4),
  warning = function(x) { w <<- conditionMessage(x)
    invokeRestart("muffleWarning") })
cat("warning:", w, "\n")
cat("identical to prior_draws():", identical(ps, prior_draws(f4)), "\n")
cat("names:", names(ps), "\n")
ps2 <- suppressWarnings(prior_samples(f4, variable = "b_x"))
cat("variable = b_x: names", names(ps2), " nrow", nrow(ps2), "\n")

cat("\n== D posterior_average(), pp_average(), model_weights()\n")
f5 <- q(y ~ 1, data = d)
r <- tryCatch(posterior_average(f1, f5, variable = c("b_x", "nope"),
                                weights = c(0.5, 0.5)),
              error = function(e) conditionMessage(e))
cat("missing var:", r, "\n")
r <- tryCatch(posterior_average(f1, f5, variable = "nope", missing = 0,
                                weights = c(0.5, 0.5)),
              error = function(e) conditionMessage(e))
cat("missing var, missing = 0:", r, "\n")
mw <- model_weights(f1, f5, weights = "waic")
cat("model_weights waic:", format(mw, digits = 6), " names", names(mw), "\n")
pa <- pp_average(f1, f5, weights = c(0.7, 0.3), seed = 3)
cat("pp_average summary: class", class(pa), " dim", dim(pa), "\n")
cat("colnames", colnames(pa), "\n")
cat("attr weights", attr(pa, "weights"), " ndraws", attr(pa, "ndraws"), "\n")
pa2 <- pp_average(f1, f5, weights = c(0.7, 0.3), seed = 3, summary = FALSE,
                  method = "posterior_epred", newdata = d[1:3, ])
cat("pp_average epred raw: dim", dim(pa2), "\n")
r <- tryCatch(pp_average(f1, f5, weights = c(0.7, 0.3), draw_ids = 1:3),
              error = function(e) conditionMessage(e))
cat("draw_ids:", r, "\n")
cat("DONE\n")
