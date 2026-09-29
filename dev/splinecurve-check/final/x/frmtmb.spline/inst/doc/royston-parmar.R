## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4.2)
has_tinyplot <- requireNamespace("tinyplot", quietly = TRUE)
has_flexsurv <- requireNamespace("flexsurv", quietly = TRUE)
knitr::opts_chunk$set(eval = has_flexsurv)
library(frmtmb)
library(frmtmb.spline)

## ----data---------------------------------------------------------------------
data(bc, package = "flexsurv")
table(bc$group, event = bc$censrec)

## ----cens---------------------------------------------------------------------
bc$censored <- 1 - bc$censrec

## ----fit----------------------------------------------------------------------
fit <- frm(bf(recyrs | cens(censored) ~ group),
           family = royston_parmar(df = 3), data = bc)
fixef(fit)

## ----hr-----------------------------------------------------------------------
exp(unlist(fixef_by_dpar(fit)$mu)[-1])

## ----identity-----------------------------------------------------------------
fs <- flexsurv::flexsurvspline(survival::Surv(recyrs, censrec) ~ group,
                               data = bc, k = 2, scale = "hazard")
kn <- fs$knots
o <- frm(bf(recyrs | cens(censored) ~ group),
         family = royston_parmar(knots = kn[2:(length(kn) - 1)],
                                 bknots = kn[c(1, length(kn))]),
         data = bc, dry_run = "objective")
cf <- fs$res[, "est"]
gam <- cf[grep("^gamma", names(cf))]
par <- c(gam[1], cf[!grepl("^gamma", names(cf))], gam[-1])
c(flexsurv = fs$loglik, frmtmb_at_flexsurv_par = -o$obj$fn(par))

## ----optima-------------------------------------------------------------------
c(flexsurv = fs$loglik, frmtmb = as.numeric(logLik(fit)))

## ----baseline-----------------------------------------------------------------
tg <- exp(seq(log(0.3), log(max(bc$recyrs)), length.out = 60))
nd <- data.frame(recyrs = tg,
                 group = factor("Good", levels = levels(bc$group)))
lp <- lapply(c("mu", "gamma1", "gamma2", "gamma3"), function(dp) {
  frm_linpred(fit, newdata = nd, type = "link", dpar = dp)
})
kn <- environment(family(fit)[["lpdf"]])$allknots
x <- log(tg)
# the Royston-Parmar basis, written out: a constant, x itself, and one
# natural cubic term per interior knot, each linear beyond the boundary
rp_basis <- function(knots, x) {
  nk <- length(knots)
  out <- list(rep(1, length(x)), x)
  for (j in seq_len(nk - 2)) {
    kj <- knots[j + 1]
    lam <- (knots[nk] - kj) / (knots[nk] - knots[1])
    out[[j + 2]] <- pmax(x - kj, 0)^3 - lam * pmax(x - knots[1], 0)^3 -
      (1 - lam) * pmax(x - knots[nk], 0)^3
  }
  out
}
bas <- rp_basis(kn, x)
logH <- Reduce(`+`, Map(function(b, g) b * as.numeric(g), bas, lp))
S <- exp(-exp(logH))
head(data.frame(t = round(tg, 3), logH = round(logH, 3),
                survival = round(S, 3)), 4)

## ----fig-surv, eval = has_tinyplot && has_flexsurv, fig.alt = "Fitted survival curves for the three prognostic groups against time in years, each falling from 1. The good-prognosis curve stays above 0.8 at eight years, the medium curve reaches about 0.6 and the poor curve about 0.3. Step functions of the same three colors, the Kaplan-Meier estimates, track each fitted curve closely."----
cols <- c("steelblue4", "goldenrod3", "firebrick")
surv_of <- function(g) {
  ndg <- data.frame(recyrs = tg,
                    group = factor(g, levels = levels(bc$group)))
  lps <- lapply(c("mu", "gamma1", "gamma2", "gamma3"), function(dp) {
    as.numeric(frm_linpred(fit, newdata = ndg, type = "link", dpar = dp))
  })
  exp(-exp(Reduce(`+`, Map(`*`, bas, lps))))
}
tinyplot::tinyplot(x = tg, y = surv_of("Good"), type = "l", lwd = 2,
                   col = cols[1], theme = "clean2", ylim = c(0, 1),
                   xlab = "years", ylab = "survival",
                   main = "Fitted survival against Kaplan-Meier")
for (i in 2:3) {
  tinyplot::tinyplot_add(x = tg, y = surv_of(levels(bc$group)[i]),
                         type = "l", lwd = 2, col = cols[i])
}
km <- survival::survfit(survival::Surv(recyrs, censrec) ~ group, data = bc)
km_s <- summary(km)
for (i in seq_along(levels(bc$group))) {
  k <- km_s$strata == levels(km_s$strata)[i]
  lines(km_s$time[k], km_s$surv[k], type = "s", col = cols[i], lty = 3)
}
legend("bottomleft", levels(bc$group), col = cols, lwd = 2, bty = "n")

## ----tvc----------------------------------------------------------------------
ftv <- frm(bf(recyrs | cens(censored) ~ group, gamma1 ~ group),
           family = royston_parmar(df = 3), data = bc)
c(proportional = as.numeric(logLik(fit)),
  time_varying = as.numeric(logLik(ftv)),
  AIC_proportional = AIC(fit), AIC_time_varying = AIC(ftv))

## ----refuse, error = TRUE-----------------------------------------------------
try({
fitted(fit)
})

## ----floored------------------------------------------------------------------
rp_floored(fit, action = "report")

## ----floored-bad--------------------------------------------------------------
set.seed(20260905)
n <- 600
bad <- data.frame(grp = factor(rep(c("A", "B"), each = n / 2)))
bad$t <- rweibull(n, shape = 1.6,
                  scale = ifelse(bad$grp == "A", 0.30, 3.0))
bad$censored <- 0L
bad$t[1] <- 50
bad$censored[1] <- 1L
bad_fit <- frm(bf(t | cens(censored) ~ grp),
               family = royston_parmar(df = 3), data = bad)
c(converged = bad_fit$opt$convergence, logLik = as.numeric(logLik(bad_fit)))
str(rp_floored(bad_fit, action = "report"))

## ----monotone-----------------------------------------------------------------
dlogH <- diff(logH) / diff(x)
c(min_slope = min(dlogH), all_increasing = all(dlogH > 0))

## ----df-----------------------------------------------------------------------
sapply(1:5, function(k) {
  AIC(frm(bf(recyrs | cens(censored) ~ group),
          family = royston_parmar(df = k), data = bc))
})

