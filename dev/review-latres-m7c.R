# Short check, item 6: M7's stopgap. Which routes lose an SE they could read at no
# cost, and what frm_curve() says. The trunc NaN of M6. Names of a two-block model.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
set.seed(1)
N <- 4000
dd <- data.frame(t = seq_len(N), g = factor(1), x = rnorm(N), f = factor(sample(letters[1:3], N, TRUE)))
dd$y <- rpois(N, exp(.5 + .3 * dd$x + as.numeric(arima.sim(list(ar = .5), N, sd = .5))))
fb <- suppressWarnings(frm(y ~ x + f + ar(t, g, cov = TRUE), family = poisson(), data = dd))
r <- msg_of(frm_curve(fb, newdata = dd[1:10, ]))
cat("frm_curve at 4000 rows:", r$err, "\n\n")
ce <- suppressMessages(conditional_effects(fb, "x"))$x
cat("conditional_effects lower__ NA:", all(is.na(ce$lower__)), " rows", nrow(ce), "\n")
fn <- suppressMessages(fitted(fb, newdata = dd[1:50, ]))
cat("fitted(newdata) Est.Error NA:", all(is.na(fn[, 2])), "\n")
op <- options(frmtmb.latres_se_rows = 1e6)
t1 <- system.time(ce2 <- conditional_effects(fb, "x")$x)[[3]]
cat("limit moved: conditional_effects", t1, "s; lower__ NA:", all(is.na(ce2$lower__)), "\n")
t2 <- system.time(fn2 <- fitted(fb, newdata = dd[1:50, ]))[[3]]
cat("limit moved: fitted(newdata)", t2, "s; Est.Error", round(fn2[1:2, 2], 4), "\n")
t3 <- system.time(lp <- frm_linpred(fb, newdata = dd[1:50, ], se.fit = TRUE))[[3]]
cat("limit moved: frm_linpred(newdata, se.fit)", t3, "s; se", round(lp$se.fit[1:2], 4), "\n")
t4 <- system.time(cv <- frm_curve(fb, newdata = dd[1:10, ]))[[3]]
cat("limit moved: frm_curve", t4, "s\n")
options(op)
cat("\n== trunc: what the node reads\n")
set.seed(167)
d <- expand.grid(t = 1:6, g = factor(1:25))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
d$cnt <- rpois(n, exp(.8 + .3 * d$x + e))
nd <- data.frame(t = 1:3, g = factor(rep("new", 3)), x = c(-1, 0, 1))
ft <- suppressWarnings(frm(cnt | trunc(ub = 100) ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
print(fitted(ft, newdata = nd, allow_new_levels = TRUE))
ft2 <- suppressWarnings(frm(cnt | trunc(lb = 0, ub = 1000) ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
print(fitted(ft2, newdata = nd, allow_new_levels = TRUE))
ft3 <- suppressWarnings(frm(cnt | trunc(lb = 0) ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
print(fitted(ft3, newdata = nd, allow_new_levels = TRUE))
rm_ <- fi("response_mean")
fam <- ft$spec$responses[[1]]$family
for (mu in c(10, 50, 100, 200, 500, 1000, 5000)) {
  r <- msg_of(rm_(fam, list(mu = mu), list(trunc_ub = 100, trunc_lb = -Inf)))
  cat("response_mean trunc mu", mu, ":", r$val %||% r$err, "\n")
}
cat("\n== two latent blocks, names\n")
suppressPackageStartupMessages(library(frmtmb.sample))
d$pos <- rgamma(n, 4, 4 / exp(.3 * d$x + e))
dm <- frm_sample(frmtmb::bf(cnt ~ x + ar(t, g, cov = TRUE)) + frmtmb::bf(pos ~ x + cosy(t, g), family = Gamma(link = "log")), family = poisson(), data = d, chains = 1, iter = 60, warmup = 30, seed = 11, refresh = 0)
vm <- variables(dm)
i <- grep("latent_state", vm)
print(vm[i[c(1, 2, 149, 150, 151, 152, length(i))]])
cat("n", length(i), " rows", n, " anyDuplicated", anyDuplicated(vm), "\n")
fm <- suppressWarnings(frm(frmtmb::bf(cnt ~ x + ar(t, g, cov = TRUE)) + frmtmb::bf(pos ~ x + cosy(t, g), family = Gamma(link = "log")), family = poisson(), data = d))
lab <- fi("brms_r_labels")(fm)
print(length(lab)); print(anyDuplicated(lab)); print(lab[c(1, 150, 151, 152)])
print(tail(names(posterior::as_draws_df(dm)), 12))
