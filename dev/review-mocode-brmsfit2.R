## Fix round 1 check, brms 2.23.0 side (rstan; 2 chains of 2000, seed 3):
## the values conditional_effects() holds two mo() variables at and the
## moderator values of a mo() variable; posterior sds of (mo(k) | g) at
## a known level, a new level and re_formula = NA. Data as in
## dev/review-mocode-ce3.R (seed 5).
##   Rscript dev/review-mocode-brmsfit2.R > dev/review-mocode-log/brmsfit2.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
f7 <- function(x) paste(format(signif(as.numeric(x), 7)), collapse = " ")
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$j <- sample(2:6, n, replace = TRUE, prob = c(1, 1, 3, 1, 1))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] +
  c(0, 0.5, 0.8, 1, 1.1)[d$j - 1L] + rnorm(8, sd = 0.4)[d$g] +
  d$k * rnorm(8, sd = 0.2)[d$g] + rnorm(n)
cat("median/mad k:", median(d$k), mad(d$k), " j:", median(d$j), mad(d$j),
    " mean k, j:", mean(d$k), mean(d$j), "\n")
f2 <- q(brm(y ~ mo(k) + mo(j) + z, data = d, chains = 2, iter = 2000,
            seed = 3, refresh = 0, backend = "rstan"))
for (eff in c("z", "k", "j", "z:k", "z:j", "k:j")) {
  ce <- q(conditional_effects(f2, effects = eff))[[1L]]
  cat(sprintf("ce(%s): k = {%s} j = {%s} z = {%s} effect2__ = {%s}\n", eff,
              f7(unique(ce$k)), f7(unique(ce$j)),
              f7(utils::head(unique(ce$z), 2)),
              paste(unique(as.character(ce$effect2__)), collapse = " ")))
  cat(sprintf("        estimate__ (first 4): %s\n", f7(head(ce$estimate__, 4))))
}
fre <- q(brm(y ~ mo(k) + (mo(k) | g), data = d, chains = 2, iter = 2000,
             seed = 3, refresh = 0, backend = "rstan"))
sdp <- function(e) apply(e, 2, sd)
cat("(mo(k) | g) g = 1, k = 1..4 sd:",
    f7(sdp(q(posterior_epred(fre, newdata = data.frame(k = 1:4, g = "1"))))),
    "\n")
cat("(mo(k) | g) re_formula = NA sd:",
    f7(sdp(q(posterior_epred(fre, newdata = data.frame(k = 1:4),
                             re_formula = NA)))), "\n")
cat("(mo(k) | g) new level sd:",
    f7(sdp(q(posterior_epred(fre, newdata = data.frame(k = 1:4, g = "new"),
                             allow_new_levels = TRUE)))), "\n")
cat("(mo(k) | g) g = 1 mean:",
    f7(colMeans(q(posterior_epred(fre, newdata = data.frame(k = 1:4,
                                                            g = "1"))))),
    "\n")
