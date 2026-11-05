# M6, measured: a difference at one x and one sdx, against brms 2.23.0 at the same
# parameter values, and against a brms fit's own draws of the same contrast. Data
# seed 7402, as dev/review-predread-regress-me.R. Also two brms readings the fix
# round needs: a missing me(gr = ) level without allow_new_levels, and a display
# grid conditioned on one level.
#   Rscript dev/predread-me-contrast.R > dev/predread-log/me-contrast.txt
source("dev/review-predread-lib.R")
set.seed(7402)
n <- 200
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.5,
                g = factor(rep(sprintf("g%02d", 1:20), each = 10)))
d$xo <- d$lat + rnorm(n, 0, 0.5)
lv <- rnorm(20)[as.integer(d$g)]
d$xg <- lv + rnorm(20, 0, 0.4)[as.integer(d$g)]
d$sg <- 0.4
d$y <- 1 + 1.2 * d$lat + 0.6 * d$z + rnorm(n, 0, 0.5)
d$y2 <- 1 + 1.1 * lv + 0.6 * d$z + rnorm(n, 0, 0.5)
fit <- quiet(frm(bf(y ~ me(xo, sx) + z), data = d, family = gaussian()))
fx <- quiet(rev_E$brms_fixed_fit(brms::bf(y ~ me(xo, sx) + z), stats::gaussian(),
                                 d, fit, ndraws = 4000))
nd1 <- data.frame(xo = c(-1, 1), sx = 0.5, z = 1)
nd0 <- data.frame(xo = c(-1, 1), sx = 0.5, z = 0)
b <- fixef(fit)
cat("coefficient of z:", round(b["z", c("Estimate", "Est.Error")], 4), "\n")
cv <- frm_curve(fit, newdata = nd1, contrast = nd0, simultaneous = FALSE)
cat("frm_curve(contrast =) estimate", round(cv$.estimate, 4), " se",
    round(cv$.se, 4), "\n")
em <- as.data.frame(pairs(emmeans::emmeans(fit, ~ z, at = list(z = c(0, 1),
                                                               xo = 1, sx = 0.5))))
cat("frmtmb emmeans pairs estimate", round(em$estimate, 4), " se", round(em$SE, 4),
    "\n")
eb <- as.data.frame(pairs(emmeans::emmeans(fx, ~ z, at = list(z = c(0, 1),
                                                              xo = 1, sx = 0.5))))
print(eb)
a <- posterior_epred(fx, newdata = nd1)
c0 <- posterior_epred(fx, newdata = nd0)
dd <- a - c0
cat("brms epred(nd1) - epred(nd0), sd over draws:", round(apply(dd, 2, sd), 4),
    " mean", round(colMeans(dd), 4), "\n")
cat("brms 95 percent of it:", round(quantile(dd[, 2], c(0.025, 0.975)), 3), "\n")
# on a fit brms sampled, one newdata carrying both rows: still one draw per row
both <- rbind(nd1, nd0)
a2 <- posterior_epred(fx, newdata = both)
cat("one call, rows 1 and 3 at one x:", round(sd(a2[, 1] - a2[, 3]), 4), "\n")
# me(gr = ): one level on both grids reads one value, in both packages
fg <- quiet(frm(bf(y2 ~ me(xg, sg, gr = g) + z), data = d, family = gaussian()))
gx <- quiet(rev_E$brms_fixed_fit(brms::bf(y2 ~ me(xg, sg, gr = g) + z),
                                 stats::gaussian(), d, fg, ndraws = 2000))
n1 <- data.frame(xg = d$xg[1], sg = 0.4, z = 1, g = d$g[1])
n0 <- transform(n1, z = 0)
cvg <- frm_curve(fg, newdata = rbind(n1, n1), contrast = rbind(n0, n0),
                 simultaneous = FALSE)
cat("gr frm_curve(contrast =) se", round(cvg$.se[1], 4), " coefficient se",
    round(fixef(fg)["z", "Est.Error"], 4), "\n")
ag <- posterior_epred(gx, newdata = rbind(n1, n0))
cat("gr brms, one newdata, sd of row 1 - row 2:", round(sd(ag[, 1] - ag[, 2]), 6),
    "\n")
ag1 <- posterior_epred(gx, newdata = n1)
ag0 <- posterior_epred(gx, newdata = n0)
cat("gr brms, two calls, sd of the difference:", round(sd(ag1 - ag0), 4), "\n")
nan <- n1
nan$g <- NA
cat("brms gr NA level, no allow_new_levels:", tryCatch(
  round(sd(posterior_epred(gx, newdata = nan)), 4),
  error = function(e) conditionMessage(e)), "\n")
ce <- tryCatch(quiet(brms::conditional_effects(gx, "xg",
                     conditions = data.frame(g = d$g[1]), resolution = 3)),
               error = function(e) conditionMessage(e))
cat("brms CE xg, conditions g =", if (is.character(ce)) ce else
  round(ce[[1L]]$estimate__, 3), "\n")
