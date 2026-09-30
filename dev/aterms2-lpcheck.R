# brms_lp_check() for rate() and subset() before the rows go into
# test-brms-likelihood.R. Seeds inside. Output:
# dev/aterms2-log-lpcheck.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/stan-cache",
           FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
options(frmtmb.brms_lp_report = TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-aterms2/tests/testthat/helper-brms.R")
run <- function(label, expr) {
  cat("\n==========", label, "\n")
  r <- tryCatch(expr, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) {
    cat(sprintf("measured_const %.6g max_grad %.3g ours %.10g\n",
                r$measured_const, r$max_grad, r$ours))
  }
  invisible(r)
}

set.seed(24)
n <- 300
dr <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
dr$y <- rpois(n, exp(0.2 + 0.5 * dr$x) * dr$time)
dr$yn <- rnbinom(n, mu = exp(0.2 + 0.5 * dr$x) * dr$time,
                 size = 2 * dr$time)
fp <- frm(y | rate(time) ~ x, data = dr, family = poisson())
run("rate poisson", brms_lp_check(brms::bf(y | rate(time) ~ x), poisson(),
                                  dr, fp))
fn <- frm(yn | rate(time) ~ x, data = dr, family = negbinomial())
run("rate negbinomial", brms_lp_check(brms::bf(yn | rate(time) ~ x),
                                      brms::negbinomial(), dr, fn))
fns <- frm(bf(yn | rate(time) ~ x, shape ~ x), data = dr,
           family = negbinomial())
run("rate negbinomial shape ~ x",
    brms_lp_check(brms::bf(yn | rate(time) ~ x, shape ~ x),
                  brms::negbinomial(), dr, fns))
fg <- frm(y | rate(time) ~ x, data = dr, family = geometric())
run("rate geometric", brms_lp_check(brms::bf(y | rate(time) ~ x),
                                    brms::geometric(), dr, fg))
fpi <- frm(y | rate(time) ~ x, data = dr, family = poisson("sqrt"))
run("rate poisson sqrt", brms_lp_check(brms::bf(y | rate(time) ~ x),
                                       poisson("sqrt"), dr, fpi))

set.seed(25)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                s1 = rep(c(TRUE, FALSE), n / 2),
                s2 = c(rep(TRUE, 150), rep(FALSE, 50)))
d$y1 <- 1 + d$x + rnorm(n)
d$y2 <- rpois(n, exp(0.5 - 0.3 * d$z))
d$y2[!d$s2] <- NA
d$z[!d$s2] <- NA
bfs <- brms::bf(y1 | subset(s1) ~ x, family = gaussian()) +
  brms::bf(y2 | subset(s2) ~ z, family = poisson()) +
  brms::set_rescor(FALSE)
fs <- frm(bf(y1 | subset(s1) ~ x) + gaussian() +
            bf(y2 | subset(s2) ~ z) + poisson(), data = d)
run("subset two responses", brms_lp_check(bfs, NULL, d, fs))

set.seed(26)
n <- 240
dm <- data.frame(g1 = sample(seq(1, n - 1, 2), n, TRUE), g2 = seq_len(n),
                 s = rep(c(TRUE, FALSE), n / 2), w = rnorm(n))
dm$x <- rnorm(n)
dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + 0.3 * dm$w + rnorm(n, sd = 0.5)
dm$x[c(3, 7, 8, 21)] <- NA
bm <- brms::bf(y ~ mi(x, idx = g1) + w) +
  brms::bf(x | mi() + index(g2) + subset(s) ~ 1) + brms::set_rescor(FALSE)
fm <- frm(bf(y ~ mi(x, idx = g1) + w) +
            bf(x | mi() + index(g2) + subset(s) ~ 1), data = dm,
          family = gaussian())
print(summary(fm))
run("subset + mi idx", brms_lp_check(bm, gaussian(), dm, fm, joint = TRUE))
cat("\nbrms prior coef names:\n")
print(brms::default_prior(bm, dm)[, c("class", "coef", "resp")])
