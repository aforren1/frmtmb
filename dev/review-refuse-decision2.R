# Item 3: what the two routes do to what a user decides, and what the observed
# rows route would break. Seed 9500. Reads no lane script.
# (a) family choice by AIC under joint Laplace and on the observed rows;
# (b) the L1 of the observed-rows route: brms's mi() program minus the Ymi
#     density is the complete-case density at any Ymi, by hand;
# (c) fitted at the missing rows: predict(cc, newdata) vs the joint fit;
# (d) a response read on some rows only: rows nothing reads could drop;
# (e) the read case's own Laplace error stays: joint fit vs a 1-D exact ML.
source("dev/review-refuse-lib.R")
gen <- function(n, what) {
  x <- rnorm(n); eta <- 0.3 + 0.5 * x
  y <- switch(what,
    exgaussian = eta + rnorm(n, 0, 0.5) + rexp(n, 1 / 2),
    student = eta + rt(n, 3),
    skew_normal = { del <- 6 / sqrt(37); eta + del * abs(rnorm(n)) + sqrt(1 - del^2) * rnorm(n) },
    gaussian = rnorm(n, eta, 1))
  data.frame(x = x, y = y)
}
fams <- list(gaussian = gaussian(), student = student(), exgaussian = exgaussian(), skew_normal = skew_normal())
R <- 30L
set.seed(9500)
for (truth in c("exgaussian", "student", "skew_normal")) for (m in c(15, 60)) {
  pick_j <- pick_c <- character(R)
  dA <- matrix(NA_real_, R, 2)
  for (r in seq_len(R)) {
    d <- gen(200, truth); d$y[sample(200, m)] <- NA
    aj <- vapply(fams, function(f) {
      a <- tr(frm(bf(y | mi() ~ x), family = f, data = d)); if (is_err(a)) NA_real_ else AIC(a)
    }, 0)
    ac <- vapply(fams, function(f) {
      a <- tr(frm(bf(y ~ x), family = f, data = d[!is.na(d$y), ])); if (is_err(a)) NA_real_ else AIC(a)
    }, 0)
    pick_j[r] <- names(fams)[which.min(aj)]
    pick_c[r] <- names(fams)[which.min(ac)]
    dA[r, ] <- c(aj[[truth]] - aj[["gaussian"]], ac[[truth]] - ac[["gaussian"]])
  }
  tb <- function(p) paste(names(table(p)), table(p), sep = "=", collapse = " ")
  pr("(a) truth %-11s m=%2d R=%d: AIC picks, joint [%s]; observed rows [%s]; true family right %d vs %d; median AIC(truth)-AIC(gaussian) joint %.1f, observed %.1f",
     truth, m, R, tb(pick_j), tb(pick_c), sum(pick_j == truth), sum(pick_c == truth), median(dA[, 1]), median(dA[, 2]))
}

# (b) by hand: brms's mi() log_prob, the Ymi density taken off, is the
# complete-case density, whatever Ymi is, on one family per chart and one
# identity family
set.seed(9501)
one <- function(lab, fb, ff, n, yfun, ld) {
  x <- rnorm(n)
  y <- yfun(x)
  d <- data.frame(x = x, y = y); rows <- seq(5, n, 10); d$y[rows] <- NA
  bb <- brms::bf(y | mi() ~ x); pri <- E$brms_flat_prior(bb, data = d, family = fb)
  code <- brms::make_stancode(bb, data = d, family = fb, prior = pri)
  sdat <- E$brms_standata(bb, data = d, family = fb, prior = pri)
  sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 0))
  cc <- frm(bf(y ~ x), family = ff, data = d[!is.na(d$y), ])
  cf <- fixef(cc)
  b0 <- cf["Intercept", 1]; b1 <- cf["x", 1]; bd <- cc$estimates$betad
  res <- vapply(1:3, function(s) {
    set.seed(s); ymi <- ld$r(length(rows))
    pars <- c(list(b = array(b1, 1L), Intercept = b0 + b1 * mean(d$x), Ymi = array(ymi)), ld$pars(bd))
    up <- rstan::unconstrain_pars(sf, pars)
    lp <- rstan::log_prob(sf, up, adjust_transform = FALSE)
    lp - sum(ld$d(ymi, b0 + b1 * x[rows], bd))
  }, 0)
  pr("(b) %-10s brms mi() lp - Ymi density, three Ymi: %.6f %.6f %.6f; frm complete-case logLik %.6f; max diff %.2e",
     lab, res[1], res[2], res[3], as.numeric(logLik(cc)), max(abs(res - as.numeric(logLik(cc)))))
}
one("Gamma", brms::brmsfamily("gamma", "log"), stats::Gamma(link = "log"), 120,
    function(x) rgamma(length(x), 2, 2 / exp(0.3 + 0.5 * x)),
    list(r = function(k) rgamma(k, 2, 1), pars = function(bd) list(shape = exp(bd[[1]])),
         d = function(y, eta, bd) dgamma(y, exp(bd[[1]]), exp(bd[[1]]) / exp(eta), log = TRUE)))
one("gaussian", brms::brmsfamily("gaussian"), stats::gaussian(), 120,
    function(x) rnorm(length(x), 0.3 + 0.5 * x, 0.7),
    list(r = function(k) rnorm(k), pars = function(bd) list(sigma = exp(bd[[1]])),
         d = function(y, eta, bd) dnorm(y, eta, exp(bd[[1]]), log = TRUE)))

# (c) fitted at the missing rows
set.seed(9502)
d <- gen(200, "exgaussian"); rows <- sort(sample(200, 60)); d$y[rows] <- NA
a <- frm(bf(y | mi() ~ x), family = exgaussian(), data = d)
b <- frm(bf(y ~ x), family = exgaussian(), data = d[-rows, ])
fa <- fitted(a)[rows, 1]; fb <- fitted(b, newdata = d[rows, ])[, 1]
pa <- predict(b, newdata = d[rows, ])
pr("(c) fitted at missing rows: joint vs observed rows, median |diff| %.4f, max %.4f; predict(observed fit, newdata) finite: %s", median(abs(fa - fb)), max(abs(fa - fb)), all(is.finite(pa[, 1])))
pr("    the chart-mode 'imputation' of the joint fit vs the family mean at the row: Gamma ratio 1 (mode on log = mean); exgaussian mode vs mean: median diff %.3f", median(frmtmb:::mi_values(a, "y")[rows] - fa))

# (d) a response read on some rows only; rows nothing reads could drop
set.seed(9503)
n <- 200
d <- gen(n, "exgaussian")
d$z <- 0.5 + 0.6 * d$y + rnorm(n, 0, 0.7)
rows <- sort(sample(n, 60))
d$y[rows] <- NA
d$z[rows[1:40]] <- NA  # 40 of the missing y are read by nothing
j <- frm(bf(y | mi() ~ x) + exgaussian() + bf(z | subset(ok) ~ mi(y)) + gaussian(), data = transform(d, ok = !is.na(z), z = ifelse(is.na(z), 0, z)))
d2 <- d[-rows[1:40], ]
j2 <- frm(bf(y | mi() ~ x) + exgaussian() + bf(z | subset(ok) ~ mi(y)) + gaussian(), data = transform(d2, ok = !is.na(z), z = ifelse(is.na(z), 0, z)))
ej <- unlist(j$estimates[c("beta", "betad")]); e2 <- unlist(j2$estimates[c("beta", "betad")])
pr("(d) read on 20 of 60 missing rows. Dropping the 40 unread: largest |change| %.4f; logLik %.3f vs %.3f (gap per dropped row %+.4f)", max(abs(ej - e2)), as.numeric(logLik(j)), as.numeric(logLik(j2)),
   (as.numeric(logLik(j)) - as.numeric(logLik(j2))) / 40)
print(round(rbind(j = ej, dropped = e2), 3))

# (e) student at 30 percent missing: one run laid open, the shape it is sold on
set.seed(9504)
d <- gen(200, "student"); rows <- sample(200, 60); d$y[rows] <- NA
a <- frm(bf(y | mi() ~ x), family = student(), data = d)
b <- frm(bf(y ~ x), family = student(), data = d[-rows, ])
nu <- function(f) 1 + exp(f$estimates$betad[["nu_(Intercept)"]])
pr("(e) true nu 3: joint nu %.2f, observed rows nu %.2f; sigma %.3f vs %.3f", nu(a), nu(b), exp(a$estimates$betad[[1]]), exp(b$estimates$betad[[1]]))
