# asym_laplace x nl on the sweep's data (dev/equiv-nl-sweep.R, seed 20261029):
# from the default starts the body stops below the linear maximum; from the
# linear fit's coefficients it reaches it. Also extension families with a
# body that has a level of its own (cox, Royston-Parmar), on dev/nlfix-extfam.R's
# data, seed 7101.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.survival)
  library(frmtmb.latent)
})
set.seed(20261029)
n <- 300
x <- round(rnorm(n), 3)
eta <- 0.3 + 0.5 * x
pos <- exp(0.4 + 0.3 * x + rnorm(n, 0, 0.4))
cnt <- rpois(n, exp(eta))
prob <- plogis(eta)
ord <- cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = 1:4)
d <- data.frame(x = x)
d$y <- eta + rnorm(n)
a <- frm(bf(y ~ x), family = asym_laplace(), data = d)
b <- suppressWarnings(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                          family = asym_laplace(), data = d))
fe <- fixef(a)[, "Estimate"]
c1 <- suppressWarnings(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                           family = asym_laplace(), data = d,
                           start = list(beta = c("b0_(Intercept)" = fe[[1]],
                                                 "b1_(Intercept)" = fe[[2]]))))
cat(sprintf("asym_laplace: linear %.4f nl default start %.4f nl from linear %.4f\n",
            logLik(a), logLik(b), logLik(c1)))

set.seed(7101)
x <- rnorm(n)
d <- data.frame(x = x, id = factor(rep(1:30, each = 10)),
                tt = rep(1:10, times = 30))
d$time <- rexp(n, exp(0.5 * x))
d$ev <- rbinom(n, 1, 0.8)
d$cc <- 1 - d$ev
d$y <- ifelse(rep(rbinom(30, 1, 0.5), each = 10) == 1, -1, 1) + 0.5 * x +
  rnorm(n, 0, 0.5)
ll <- function(expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat("  ERROR:", substr(conditionMessage(r), 1, 160), "\n")
    return(NA_real_)
  }
  as.numeric(logLik(r))
}
cat("cox linear", ll(frm(bf(time | cens(cc) ~ x), family = cox(), data = d)),
    "nl with b0", ll(frm(bf(time | cens(cc) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                         family = cox(), data = d)), "\n")
cat("rp linear", ll(frm(bf(time | cens(cc) ~ x), family = royston_parmar(), data = d)),
    "nl with b0", ll(frm(bf(time | cens(cc) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
                         family = royston_parmar(), data = d)), "\n")
cat("hmm nl:\n")
ll(frm(bf(y | index(id, tt) ~ b1 * x, b1 ~ 1, nl = TRUE),
       family = frmtmb.latent::hmm(2, gaussian()), data = d))
