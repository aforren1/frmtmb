# brms 2.23.0 under sample_prior = "only": which terms of the model
# block stay outside `if (!prior_only)` for me(), mi() and a
# covariance-form ar(). No compile.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- data.frame(x = rnorm(60), sdx = 0.3, g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), t = rep(1:12, 5))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$xm <- d$x
d$xm[1:5] <- NA
model_block <- function(sc) {
  sc <- strsplit(sc, "\n")[[1]]
  i <- grep("^model", sc)
  j <- grep("^generated quantities", sc)
  cat(sc[i:(j - 1)], sep = "\n")
}
pr <- c(prior(normal(0, 1), class = b), prior(normal(0, 1), class = meanme),
        prior(normal(0, 1), class = sdme))
cat("===== me =====\n")
model_block(stancode(y ~ me(x, sdx), data = d, prior = pr,
                     sample_prior = "only"))
cat("===== mi =====\n")
model_block(stancode(bf(y ~ mi(xm)) + bf(xm | mi() ~ 1) + set_rescor(FALSE),
                     data = d,
                     prior = c(prior(normal(0, 1), class = b, resp = y),
                               prior(normal(0, 1), class = Intercept,
                                     resp = xm)),
                     sample_prior = "only"))
cat("===== ar cov =====\n")
model_block(stancode(y ~ 1 + ar(time = t, gr = h, cov = TRUE), data = d,
                     prior = prior(normal(0, 0.5), class = ar),
                     sample_prior = "only"))
