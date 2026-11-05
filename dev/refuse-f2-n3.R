# N3: the review's table (dev/review-refuse-f1-vm.R, seed 10600, one prior).
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
quiet <- function(e) suppressMessages(suppressWarnings(e))
set.seed(10600)
n <- 200; x <- rnorm(n)
d <- data.frame(x = x, y = 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2))
rows <- sample(n, 60); d$y[rows] <- NA; dc <- d[-rows, ]
pf <- set_prior("normal(0, 5)", class = "b") + set_prior("normal(0, 5)", class = "Intercept") +
  set_prior("exponential(1)", class = "sigma") + set_prior("exponential(0.5)", class = "beta")
go <- function(x, ...) quiet(frm_sample(x, prior = pf, chains = 4, iter = 1500, seed = 1,
                                        refresh = 0, cores = 1, ...))
s1 <- go(bf(y | mi() ~ x), family = exgaussian(), data = d, laplace = TRUE)
s2 <- go(bf(y ~ x), family = exgaussian(), data = dc, laplace = FALSE)
s3 <- go(frm(bf(y | mi() ~ x), family = exgaussian(), data = d), laplace = TRUE)
sm <- function(lab, s) {
  dd <- as.data.frame(posterior::as_draws_df(s))
  cat(sprintf("%-38s sigma %.3f (sd %.3f); beta %.3f (sd %.3f); b_x %.3f\n", lab,
              mean(dd$sigma), sd(dd$sigma), mean(dd$beta), sd(dd$beta), mean(dd$b_x)))
}
sm("frm_sample(formula, laplace = TRUE)", s1)
sm("frm_sample(fit, laplace = TRUE)", s3)
sm("observed rows, full sampling", s2)
