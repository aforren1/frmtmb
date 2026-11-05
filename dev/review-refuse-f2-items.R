# Fix round 2, items 4 and 5. Seed 11800.
source("dev/review-refuse-lib.R")
suppressPackageStartupMessages(library(frmtmb.sample))
quiet <- function(e) suppressMessages(suppressWarnings(e))
Gam <- stats::Gamma(link = "log"); Gb <- brms::brmsfamily("gamma", "log")

# 4a. the row count: rows with no density anywhere
set.seed(11800)
n <- 120; x <- rnorm(n); y0 <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * x))
d <- data.frame(x = x, y = y0, id = seq_len(n)); d$z <- 0.5 + 0.6 * y0 + rnorm(n, 0, 0.7); d$w <- rnorm(n)
M <- seq(4, n, by = 12); d$y[M] <- NA
a <- frm(bf(y | mi() ~ x) + Gam + bf(w ~ x) + gaussian(), data = d)
pr("y unread beside an observed w: attr %s (w keeps every row, so 120); frame rows %d", attr(logLik(a), "nobs"), nobs(a))
dz <- d; dz$z[M] <- NA; dz$w[M[1:5]] <- NA
a <- frm(bf(y | mi() ~ x) + Gam + bf(z | mi() ~ mi(y)) + gaussian(), data = dz)
b <- frm(bf(y ~ x) + Gam + bf(z ~ y) + gaussian(), data = dz[-M, ])
pr("chain, y and z gone on 10: attr %s, BIC %.4f; hand deletion attr %s, BIC %.4f", attr(logLik(a), "nobs"), BIC(a), attr(logLik(b), "nobs"), BIC(b))
a <- frm(bf(y | mi() ~ x) + Gam + bf(w | mi() ~ x) + gaussian(), data = dz)
pr("y gone on 10, w gone on 5 of them: attr %s (want 115)", attr(logLik(a), "nobs"))
dd <- d; dd$w[seq(5, n, by = 12)] <- NA
a <- frm(bf(y | mi() ~ x) + Gam + bf(w | mi() ~ x) + gaussian(), data = dd)
pr("y gone on 10, w gone on 10 others: attr %s (want 120: every row has a density)", attr(logLik(a), "nobs"))
ds <- d; ds$ok <- rep(c(TRUE, TRUE, FALSE), 40); ds$z2 <- ds$z
a <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z2 | subset(ok) ~ x) + gaussian(), data = ds))
pr("with a subset on the other response: %s; attr %s", msg(a, 100), if (is_err(a)) "-" else attr(logLik(a), "nobs"))
if (!is_err(a)) {
  pr("  row count by hand: rows with y observed or z2 in the subset: %d", sum(!is.na(ds$y) | ds$ok))
}

# 4b. can a dropped row feed an objective? Move mi_values and the placeholder
# and read the objective at one point
ori_p <- get("mi_placeholder", asNamespace("frmtmb"))
mv <- function(lab, expr) {
  a <- tr(expr())
  assignInNamespace("mi_placeholder", function(obs) stats::quantile(obs, 0.93, names = FALSE), "frmtmb")
  b <- tr(expr())
  assignInNamespace("mi_placeholder", ori_p, "frmtmb")
  if (is_err(a) || is_err(b)) { pr("  %-36s %s | %s", lab, msg(a, 80), msg(b, 80)); return(invisible()) }
  pr("  %-36s latent %d; logLik %.5f; placeholder moves it %.1e", lab, length(a$estimates$miss), as.numeric(logLik(a)), as.numeric(logLik(a)) - as.numeric(logLik(b)))
}
pr("placeholder against an objective:")
mv("chain", function() frm(bf(y | mi() ~ x) + Gam + bf(z | mi() ~ mi(y)) + gaussian(), data = dz))
mv("reader weight 0", function() frm(bf(y | mi() ~ x) + Gam + bf(z | weights(wz) ~ mi(y)) + gaussian(), data = transform(d, wz = ifelse(is.na(y), 0, 1))))
cnt0 <- rpois(n, 2)
mv("weight 0 under a log link", function() frm(bf(y | mi() ~ x) + Gam + bf(cnt | weights(wz) ~ mi(y)) + poisson(), data = transform(d, wz = ifelse(is.na(y), 0, 1), cnt = cnt0)))
mv("weight 0, dpar of the reader", function() frm(bf(y | mi() ~ x) + Gam + bf(z | weights(wz) ~ x, sigma ~ mi(y)) + gaussian(), data = transform(d, wz = ifelse(is.na(y), 0, 1))))
mv("weight 0 in a mixture reader", function() frm(bf(y | mi() ~ x) + Gam + bf(z | weights(wz) ~ 1, mu1 ~ mi(y), mu2 ~ 1) + mixture(gaussian(), gaussian(), order = "none"), data = transform(d, wz = ifelse(is.na(y), 0, 1))))
mv("idx reader, subset", function() {
  dq <- d; dq$ok <- TRUE; dq$ok[M[1:6]] <- FALSE; dq$z2 <- ifelse(dq$ok, dq$z, 0)
  frm(bf(y | mi() + index(id) ~ x) + Gam + bf(z2 | subset(ok) ~ mi(y, idx = id)) + gaussian(), data = dq)
})
mv("reader with a grouped RE", function() frm(bf(y | mi() ~ x) + Gam + bf(z | weights(wz) ~ mi(y) + (1 | g)) + gaussian(), data = transform(d, wz = ifelse(is.na(y), 0, 1), g = factor(rep(1:12, 10)))))
mv("reader with RE slope on another column", function() frm(bf(y | mi() ~ x) + Gam + bf(z | mi() ~ mi(y) + (0 + x | g)) + gaussian(), data = transform(dz, g = factor(rep(1:12, 10)))))
# a reader weight of 0 under autocorrelation, where the lag reads the row
da <- data.frame(g = rep(1:12, each = 10), t0 = rep(1:10, 12), x = x, y = y0, zz = arima.sim(list(ar = 0.5), n) + 0.5 * y0)
da$y[M] <- NA
mv("reader with ar(cov = FALSE)", function() frm(bf(y | mi() ~ x) + Gam + bf(zz ~ mi(y) + ar(time = t0, gr = g, p = 1)) + gaussian(), data = da))
da2 <- da; da2$zz[M] <- NA
mv("reader mi + ar(cov = FALSE)", function() frm(bf(y | mi() ~ x) + Gam + bf(zz | mi() ~ mi(y) + ar(time = t0, gr = g, p = 1)) + gaussian(), data = da2))
a <- tr(frm(bf(y | mi() ~ x) + Gam + bf(zz | mi() ~ mi(y) + ar(time = t0, gr = g, p = 1)) + gaussian(), data = da2))
if (!is_err(a)) pr("  ar reader: latent %d (the lag reads z's row, and z's mean reads y, so 20 stay)", length(a$estimates$miss))
da3 <- da; da3$s <- da3$t0 != 5
a <- tr(frm(bf(y | mi() ~ x) + Gam + bf(zz | subset(s) ~ mi(y, idx = t0) + x) + gaussian(), data = da3))
pr("  subset with mi idx nonunique: %s", msg(a, 100))

# 4c. a plug-in: fitted(z) at a row that left, with a nonlinear reader, against
# brms's posterior mean of epred at that row
set.seed(11801)
n <- 150; x <- rnorm(n)
y0 <- rgamma(n, 1.2, 1.2 / exp(0.2 + 0.4 * x))
d <- data.frame(x = x, y = y0)
d$cnt <- rpois(n, exp(0.2 + 0.5 * y0)); d$z <- 0.5 + 0.6 * y0 + rnorm(n, 0, 0.5)
M <- seq(4, n, by = 10); d$y[M] <- NA; d$w <- ifelse(is.na(d$y), 0, 1)
d$z[M] <- NA
fa <- frm(bf(y | mi() ~ x) + Gam + bf(cnt | weights(w) ~ mi(y)) + poisson(), data = d)
pr("poisson reader, weight 0 on the missing rows: latent %d", length(fa$estimates$miss))
fv <- fitted(fa, resp = "cnt")[M, 1]
bb <- brms::bf(y | mi() ~ x, family = Gb) + brms::bf(cnt | weights(w) ~ mi(y), family = poisson()) + brms::set_rescor(FALSE)
code <- brms::make_stancode(bb, data = d)
sdat <- E$brms_standata(bb, data = d)
sf <- quiet(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 4, iter = 2000, seed = 11, refresh = 0, cores = 1))
e <- rstan::extract(sf)
# brms's epred at the row: exp(Intercept + b * Ymi), by its own centering
Xc <- mean(c(d$y[!is.na(d$y)]))
yl <- matrix(rep(sdat$Y_y, each = nrow(e$Ymi_y)), nrow(e$Ymi_y))
yl[, sdat$Jmi_y] <- e$Ymi_y
mu <- exp(as.numeric(e$Intercept_cnt) + as.numeric(e$bsp_cnt) * yl)
bm <- colMeans(mu)[M]
pr("poisson reader: fitted at left rows: frm plug-in median %.3f; brms posterior mean of epred median %.3f; relative gap median %.3f, max %.3f; on observed rows median %.3f",
   median(fv), median(bm), median(abs(fv - bm) / bm), max(abs(fv - bm) / bm), median(abs(fitted(fa, resp = "cnt")[-M, 1] - colMeans(mu)[-M]) / colMeans(mu)[-M]))
# at frm's own parameters: E[exp(b0 + b1 Y)] under Y ~ Gamma, by integration
cf <- fixef(fa); b0 <- cf["cnt_Intercept", 1]; b1 <- cf["cnt_miy", 1]
shape <- exp(fa$estimates$betad[[1]])
muy <- exp(cf["y_Intercept", 1] + cf["y_x", 1] * d$x[M])
exact_m <- vapply(muy, function(m) integrate(function(v) exp(b0 + b1 * v) * dgamma(v, shape, shape / m), 0, 80 * m)$value, 0)
pr("  at frm's parameters, E[exp(b0 + b1 Y)] by integrate: median %.3f; plug-in %.3f; ratio max %.3f", median(exact_m), median(fv), max(exact_m / fv))
# the linear reader: the plug-in is the mean
fl <- frm(bf(y | mi() ~ x) + Gam + bf(z | mi() ~ mi(y)) + gaussian(), data = d)
cf <- fixef(fl)
fz <- fitted(fl, resp = "z")[M, 1]; fy <- fitted(fl, resp = "y")[M, 1]
pr("linear reader: fitted(z) - (b0 + b1 fitted(y)) at left rows, max %.2e", max(abs(fz - cf["z_Intercept", 1] - cf["z_miy", 1] * fy)))

# 5. N7 and N8 against brms
set.seed(11802)
n <- 120; lev <- paste0("L", 1:8)
dm <- data.frame(g1 = sample(lev, n, TRUE), g2 = sample(lev, n, TRUE), x = rnorm(n), y = rnorm(n), w1 = runif(n), w2 = runif(n), v1 = runif(n), v2 = runif(n))
dm$W <- cbind(dm$w1, dm$w2); dm$g <- factor(dm$g1); dm$h <- factor(dm$g2)
fr <- function(lab, txt) {
  a <- tr(frm(eval(parse(text = txt), envir = asNamespace("frmtmb")), data = dm, dry_run = "frame"))
  b <- tr(brms::make_standata(eval(parse(text = txt), envir = asNamespace("brms")), data = dm))
  pr("%-38s frm %-58s brms %s", lab, msg(a, 56), msg(b, 60))
}
fr("N7 W beside cbind", "bf(y ~ x + (1 | p | mm(g1, g2, weights = W)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
fr("N7 swapped columns", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w2, w1))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
fr("N7 proportional", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(2 * w1, 2 * w2))))")
fr("N7 scale F same weights", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2), scale = FALSE)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2), scale = FALSE)))")
fr("N7 scale differs", "bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2), scale = FALSE)), sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2))))")
fr("N8 id, other group", "bf(y ~ x + (1 | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(h, id = 'p')))")
fr("N8 |p|, other group", "bf(y ~ x + (1 | p | gr(g, cor = FALSE)), sigma ~ (1 | p | h))")
fr("N8 id F, id F other group", "bf(y ~ x + (1 | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(h, id = 'p', cor = FALSE)))")
fr("N8 |p| F, |p| F other group", "bf(y ~ x + (1 | p | gr(g, cor = FALSE)), sigma ~ (1 | p | gr(h, cor = FALSE)))")
fr("N8 id mm, other group", "bf(y ~ x + (1 | mm(g1, g2, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'p')))")
fr("N4 mm F+T", "bf(y ~ x + (1 | mm(g1, g2, id = 'p', cor = FALSE)), sigma ~ (1 | mm(g1, g2, id = 'p')))")
fr("N4 mm T+F", "bf(y ~ x + (1 | mm(g1, g2, id = 'p')), sigma ~ (1 | mm(g1, g2, id = 'p', cor = FALSE)))")
fr("N4 mm |p| T+F", "bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ (1 | p | mm(g1, g2, cor = FALSE)))")
fr("N4 mm all F: fits", "bf(y ~ x + (1 | mm(g1, g2, id = 'p', cor = FALSE)), sigma ~ (1 | mm(g1, g2, id = 'p', cor = FALSE)))")
# the |p| spelling of cor = FALSE drops the key in brms too? read cor_1
for (txt in c("bf(y ~ x + (1 | p | gr(g, cor = FALSE)), sigma ~ (1 | p | g))", "bf(y ~ x + (1 | p | g), sigma ~ (1 | p | gr(g, cor = FALSE)))")) {
  co <- brms::make_stancode(eval(parse(text = txt), envir = asNamespace("brms")), data = dm)
  a <- tr(frm(eval(parse(text = txt), envir = asNamespace("frmtmb")), data = dm))
  pr("%-64s brms cor_1 %s, z blocks %d; frm %s %s", txt, grepl("cor_1", co), length(gregexpr("z_[0-9]+;", co)[[1]]), msg(a, 70), if (is_err(a)) "" else paste(vapply(a$frame$re_blocks, function(k) paste0(k$covstruct, k$dim), ""), collapse = ","))
}
# importance and vcov(cluster =) on all-unread
set.seed(11803)
n <- 150; x <- rnorm(n); g <- factor(rep(1:15, 10))
d <- data.frame(x = x, g = g); d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * x + rnorm(15, 0, 0.4)[g])); M <- seq(4, n, by = 12); d$y[M] <- NA
i1 <- tr(frm(bf(y | mi() ~ x + (1 | g)), family = Gam, data = d, importance = 100)); i2 <- tr(frm(bf(y ~ x + (1 | g)), family = Gam, data = d[-M, ], importance = 100))
pr("importance = 100: %s | %s; logLik gap %s; coef gap %s", msg(i1, 80), msg(i2, 80), if (is_err(i1) || is_err(i2)) "-" else format(as.numeric(logLik(i1)) - as.numeric(logLik(i2)), digits = 3),
   if (is_err(i1) || is_err(i2)) "-" else format(max(abs(i1$estimates$beta - i2$estimates$beta)), digits = 3))
a <- frm(bf(y | mi() ~ x), family = Gam, data = d); b <- frm(bf(y ~ x), family = Gam, data = d[-M, ])
v1 <- tr(vcov(a, cluster = ~ g)); v2 <- tr(vcov(b, cluster = ~ g))
pr("vcov(cluster = ~ g): %s | %s; gap %s", msg(v1, 80), msg(v2, 80), if (is_err(v1) || is_err(v2)) "-" else format(max(abs(as.matrix(v1) - as.matrix(v2))), digits = 3))
v3 <- tr(vcov(a, cluster = rep(1:15, 10))); pr("vcov(cluster = vector of all rows): %s", msg(v3, 130))
d2 <- d; d2$z <- rnorm(n)
ra <- frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y)) + gaussian(), data = d2)
v4 <- tr(vcov(ra, cluster = ~ g)); pr("vcov(cluster) with read rows refuses: %s", msg(v4, 130))
sv <- tr(sandwich::sandwich); s1 <- tr(vcov(a, type = "HC0")); pr("HC0: %s", msg(s1, 90))
