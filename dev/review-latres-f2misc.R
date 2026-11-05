# Short check 2: the post-fit refusal's cost and text; the translator with L read
# directly; R7 beyond latent residuals; frm_lp_basis() on the fitted rows.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
options(warn = -1)
cat("== the refusal: time spent first, and the text\n")
for (n in c(150, 2000, 20000)) {
  set.seed(1)
  G <- max(1, n / 6)
  d <- data.frame(t = rep(1:6, length.out = n), g = factor(rep(seq_len(ceiling(G)), each = 6, length.out = n)))
  d$x <- rnorm(n)
  e <- rnorm(n, 0, .5)
  d$y <- rgamma(n, 2, 2 * (1.2 + .1 * d$x + e * .5 * (abs(e) < 1)))
  t0 <- proc.time()[[3]]
  r <- msg_of(frm(y ~ x + ar(t, g, cov = TRUE), family = Gamma(), data = d))
  cat(sprintf("n %d: %.1f s, %s\n", n, proc.time()[[3]] - t0, if (is.null(r$err)) "fits" else substr(r$err, 1, 430)))
  t0 <- proc.time()[[3]]
  r2 <- msg_of(frm(y ~ x + ar(t, g, cov = TRUE), family = Gamma(), data = d, dry_run = "frame"))
  cat(sprintf("   dry_run frame %.2f s: %s\n", proc.time()[[3]] - t0, if (is.null(r2$err)) "passes" else substr(r2$err, 1, 100)))
}
cat("\n== other routes through a fit that would refuse\n")
set.seed(1)
d <- expand.grid(t = 1:6, g = factor(1:25))
d$x <- rnorm(nrow(d))
d$z <- rnorm(nrow(d))
d$y <- rpois(nrow(d), 2 + .5 * d$x + rep(rnorm(25, 0, .5), each = 6))
fk <- msg_of(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson("identity"), data = d))
cat("poisson identity mean 2:", substr(fk$err %||% "fits", 1, 120), "\n")
set.seed(2)
d2 <- d
d2$y <- rpois(nrow(d), 8 + .5 * d$x + as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = 1.2))))
fk2 <- msg_of(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson("identity"), data = d2))
cat("poisson identity mean 8:", substr(fk2$err %||% "fits", 1, 120), "\n")
if (!is.null(fk2$val)) {
  f <- fk2$val
  u <- msg_of(update(f, . ~ . + z))
  cat("update:", substr(u$err %||% "fits", 1, 100), "\n")
  dd <- msg_of(drop1(f))
  cat("drop1:", substr(dd$err %||% "fits", 1, 100), "\n")
  yy <- simulate(f, nsim = 1, seed = 1)[[1]]
  rr <- msg_of(refit(f, newresp = yy))
  cat("refit:", substr(rr$err %||% "fits", 1, 100), "\n")
}

cat("\n== R7 beyond latent residuals\n")
set.seed(3)
n <- 300
dd <- data.frame(x = rnorm(n), g = factor(rep(1:30, each = 10)))
dd$pos <- rlnorm(n, .5 + .3 * dd$x, .8)
dd$zi <- rpois(n, exp(.5 + .3 * dd$x)) * rbinom(n, 1, .7)
dd$hg <- rgamma(n, 3, 3 / exp(.3 * dd$x)) * rbinom(n, 1, .7)
nd <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(dd$g)))
for (cc in list(list("lognormal", pos ~ x, lognormal()), list("zero_inflated_poisson", zi ~ x, zero_inflated_poisson()), list("hurdle_gamma", hg ~ x, hurdle_gamma()), list("poisson + (1|g)", zi ~ x + (1 | g), poisson()), list("lognormal + (1|g)", pos ~ x + (1 | g), lognormal()))) {
  f <- suppressWarnings(frm(cc[[2]], family = cc[[3]], data = dd))
  a <- msg_of(marginaleffects::predictions(f, newdata = nd))
  b <- fitted(f, newdata = nd)
  cat(sprintf("%-24s predictions %s SE %s | fitted %s SE %s  ratio %s\n", cc[[1]], paste(round(a$val$estimate, 3), collapse = " "), paste(round(a$val$std.error, 4), collapse = " "), paste(round(b[, 1], 3), collapse = " "), paste(round(b[, 2], 4), collapse = " "), paste(round(a$val$std.error / b[, 2], 3), collapse = " ")))
  cat("   get_coef names:", paste(names(marginaleffects::get_coef(f)), collapse = " "), "\n")
}
set.seed(167)
dl <- expand.grid(t = 1:6, g = factor(1:25)); dl <- dl[order(dl$g, dl$t), ]; dl$x <- rnorm(nrow(dl))
e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
dl$cnt <- rpois(nrow(dl), exp(.8 + .3 * dl$x + e))
fl <- frm(cnt ~ x + ar(t, g, cov = TRUE), family = poisson(), data = dl)
nl <- data.frame(t = 1:3, g = factor(rep("new", 3)), x = c(-1, 0, 1))
a <- marginaleffects::predictions(fl, newdata = nl)
b <- fitted(fl, newdata = nl, allow_new_levels = TRUE)
cat("latent poisson: predictions SE", round(a$std.error, 4), " fitted SE", round(b[, 2], 4), "\n")

cat("\n== frm_lp_basis on fitted rows at 4000 rows: time and memory\n")
set.seed(1)
N <- 4000
dx <- data.frame(t = seq_len(N), g = factor(1), x = rnorm(N))
dx$y <- rpois(N, exp(.5 + .3 * dx$x + as.numeric(arima.sim(list(ar = .5), N, sd = .5))))
fb <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = dx))
t1 <- system.time(z <- msg_of(frm_lp_basis(fb, newdata = dx[1:10, ])))[[3]]
cat("frm_lp_basis newdata 10 rows", t1, "s", substr(z$err %||% "ok", 1, 100), "\n")
t2 <- system.time(z <- msg_of(frm_lp_basis(fb)))[[3]]
cat("frm_lp_basis fitted rows", t2, "s", substr(z$err %||% "ok", 1, 100), "\n")
cat("help mentions cost:", any(grepl("dense|cost|expensive", readLines("man/frm_lp_basis.Rd"), ignore.case = TRUE)), "\n")
print(grep("dense|cost|expensive", readLines("man/frm_lp_basis.Rd"), ignore.case = TRUE, value = TRUE))
