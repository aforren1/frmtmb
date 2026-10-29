# Lane meme: where brms 2.23.0 writes the me() prior statements among
# the others under sample_prior = "yes", and brms's prior-draw names
# after renaming (one small fit each). Output: dev/meme-out/brms-order.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
set.seed(11)
n <- 40
d <- data.frame(x = rnorm(n), z = rnorm(n), sx = 0.3,
                sz = runif(n, 0.2, 0.4), g = factor(rep(1:8, each = 5)),
                h = factor(rep(1:5, 8)))
d$xg <- rnorm(8)[d$g]
d$sg <- 0.25
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$y2 <- d$y + rnorm(n)
gq <- function(sc) {
  l <- strsplit(as.character(sc), "\n")[[1]]
  cat(grep("^  real prior_|^  .*prior_.*_rng", l, value = TRUE), sep = "\n")
}
pr <- c(prior(normal(0, 1), class = b),
        prior(normal(0, 1), class = sd),
        prior(normal(0, 2), class = meanme),
        prior(normal(1, 0.5), class = meanme, coef = mez),
        prior(exponential(1), class = sdme),
        prior(lkj(2), class = corme))
cat("==== F2: me + me + me(gr) + (1 | h) + sigma ~ z\n")
f <- bf(y ~ me(x, sx) + me(z, sz) + me(xg, sg, gr = g) + (1 | h),
        sigma ~ z)
gq(stancode(f, data = d, prior = pr, sample_prior = "yes"))
cat("==== MV with rescor\n")
fm <- bf(mvbind(y, y2) ~ me(x, sx) + (1 | h)) + set_rescor(TRUE)
gq(stancode(fm, data = d, prior = c(
  prior(normal(0, 1), class = b, resp = y),
  prior(normal(0, 1), class = b, resp = y2),
  prior(normal(0, 2), class = meanme),
  prior(exponential(1), class = sdme)), sample_prior = "yes"))
cat("==== 'only' refusal order, b proper, meanme flat\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                       sample_prior = "only",
                       prior = c(prior(normal(0, 1), class = b))),
              error = function(e) conditionMessage(e))
print(r)
cat("==== 'only' refusal, b and meanme proper, sdme flat\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                       sample_prior = "only",
                       prior = c(prior(normal(0, 1), class = b),
                                 prior(normal(0, 1), class = meanme))),
              error = function(e) conditionMessage(e))
print(r)
cat("==== 'only' refusal, meanme only for mez\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                       sample_prior = "only",
                       prior = c(prior(normal(0, 1), class = b),
                                 prior(normal(0, 1), class = meanme,
                                       coef = mez),
                                 prior(exponential(1), class = sdme))),
              error = function(e) conditionMessage(e))
cat(if (is.character(r) && length(r) == 1) r else
  grep("meanme", strsplit(as.character(r), "\n")[[1]], value = TRUE),
  sep = "\n")
cat("==== 'only' refusal, corme flat\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                       sample_prior = "only",
                       prior = c(prior(normal(0, 1), class = b),
                                 prior(normal(0, 1), class = meanme),
                                 prior(exponential(1), class = sdme),
                                 set_prior("", class = "corme"))),
              error = function(e) conditionMessage(e))
cat(if (is.character(r) && length(r) == 1) r else
  grep("Lme", strsplit(as.character(r), "\n")[[1]], value = TRUE),
  sep = "\n")

cat("==== fit: prior draw names after renaming\n")
f1 <- suppressWarnings(suppressMessages(brm(
  bf(y ~ me(x, sx) + me(z, sz) + me(xg, sg, gr = g)), data = d,
  prior = pr[-2], sample_prior = "yes", chains = 1, iter = 200, seed = 1,
  refresh = 0)))
cat("variables:", variables(f1), "\n", fill = 78)
cat("prior_draws:", names(prior_draws(f1)), "\n")
cat("prior_summary:\n")
print(as.data.frame(prior_summary(f1))[, c("prior", "class", "coef",
                                           "group", "lb", "ub",
                                           "source")])
cat("==== fit: class-wide only\n")
f2 <- suppressWarnings(suppressMessages(brm(
  bf(y ~ me(x, sx) + me(z, sz)), data = d,
  prior = c(prior(normal(0, 2), class = meanme),
            prior(exponential(1), class = sdme)),
  sample_prior = "yes", chains = 1, iter = 200, seed = 1, refresh = 0)))
cat("prior_draws:", names(prior_draws(f2)), "\n")
print(round(sapply(prior_draws(f2), function(v) c(mean(v), sd(v))), 2))
cat("hypothesis meanme_mex = 0:\n")
print(hypothesis(f2, "meanme_mex = 0", class = NULL)$hypothesis)
