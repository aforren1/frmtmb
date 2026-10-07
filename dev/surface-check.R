# Lane surface: read the new outputs on the lane build and check them
# against an independent computation.
#
#   Rscript dev/surface-check.R > dev/surface-out/check.txt
source("dev/surface-env.R")
surface_env("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb from", find.package("frmtmb"), "\n")
options(mc.cores = 1)

## frm_multiple(): summary, and the pooled curve against a hand pooling
data("nhanes", package = "mice")
imp <- mice::mice(nhanes, m = 5, print = FALSE, seed = 1)
fm <- frm_multiple(bmi ~ age * chl, data = imp)
print(summary(fm))
ce <- conditional_effects(fm, "age:chl")
d <- ce[[1]]
# by hand: each imputation's curve on the first one's grid, through
# frm_linpred() (identity link, so the link-scale pair is the curve's)
nd <- d[c("age", "chl")]
P <- sapply(fm$fits, function(f) frm_linpred(f, newdata = nd,
                                             type = "link", se.fit = TRUE)$fit)
S <- sapply(fm$fits, function(f) frm_linpred(f, newdata = nd,
                                             type = "link", se.fit = TRUE)$se.fit)
pl <- frmtmb:::rubin_pool(P, S^2, df.residual(fm$fits[[1]]))
cat("pooled CE vs hand pooling: max |estimate diff|",
    format(max(abs(d$estimate__ - pl$estimate))), " max |se diff|",
    format(max(abs(d$se__ - pl$se))), " rows", nrow(d), "\n")
cat("pooled CE se / first imputation's se: median",
    format(median(d$se__ / conditional_effects(fm$fits[[1]],
                                                "age:chl")[[1]]$se__)), "\n")
# fixef() of the pooled fit against mice's own pooling of lm()
ref <- summary(mice::pool(with(imp, lm(bmi ~ age * chl))))
fx <- fixef(fm)
cat("fixef(fm) Estimate vs mice::pool(lm): max |diff|",
    format(max(abs(fx[, "Estimate"] - ref$estimate))), "\n")

## ordinal-family pooled curve: the category probabilities
set.seed(7)
n <- 150
ord <- lapply(1:3, function(i) {
  x <- rnorm(n)
  y <- cut(x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
  data.frame(x = x, y = y)
})
fo <- frm_multiple(y ~ x, family = cumulative(), data = ord)
print(fixef(fo))
co <- conditional_effects(fo, "x")
cat("ordinal pooled CE: rows", nrow(co[[1]]), "; probabilities in [0,1]:",
    all(co[[1]]$lower__ >= 0 & co[[1]]$upper__ <= 1), "\n")

## update() and the stale prior: the message, and the fit against a
## direct fit with the remaining priors
kidney <- brms::kidney
pr <- c(set_prior("normal(0,5)", class = "b"),
        set_prior("cauchy(0,2)", class = "sd"),
        set_prior("lkj(2)", class = "cor"))
f1 <- frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal(), prior = pr)
msg <- character()
u <- withCallingHandlers(
  update(f1, formula. = ~ . - (1 + age | patient) + (1 | patient)),
  message = function(m) {
    msg <<- c(msg, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
cat("update() message:", msg, "\n")
dd <- frm(time | cens(censored) ~ age * sex + disease + (1 | patient),
          data = kidney, family = lognormal(),
          prior = c(set_prior("normal(0,5)", class = "b"),
                    set_prior("cauchy(0,2)", class = "sd")))
cat("update() vs direct fit with the two kept priors: logLik diff",
    format(as.numeric(logLik(u)) - as.numeric(logLik(dd))), "\n")
# a prior the update itself passes is still checked
r <- tryCatch(update(f1, formula. = ~ . - (1 + age | patient) + (1 | patient),
                     prior = pr), error = conditionMessage)
cat("update(prior = <with cor>):", substr(r, 1, 90), "\n")

## pp_mixture() on a fit against the posterior of the same model
set.seed(4)
dm <- data.frame(y = c(rnorm(150, -1.5), rnorm(150, 1.5)))
fmx <- frm(y ~ 1, family = mixture(gaussian(), gaussian()), data = dm)
pf <- pp_mixture(fmx)
ds <- suppressWarnings(frm_sample(fmx, chains = 2, iter = 1500,
                                  refresh = 0, seed = 11))
pd <- pp_mixture(ds)
r <- pf[, "Est.Error", 1] / pd[, "Est.Error", 1]
cat("pp_mixture Est.Error, fit (delta) / draws (posterior SD), K = 1,",
    "rows with draws SD > 0.01:\n")
print(summary(r[pd[, "Est.Error", 1] > 0.01]))
cat("pp_mixture Estimate, fit vs draws: max |diff|",
    format(max(abs(pf[, "Estimate", 1] - pd[, "Estimate", 1]))), "\n")
cat("dimnames agree:", identical(dimnames(pf)[[3]], dimnames(pd)[[3]]), "\n")

## add_criterion() on draws, brms's use_stored
ds2 <- add_criterion(ds, c("loo", "waic"))
cat("add_criterion stored:", names(ds2$criteria),
    "| loo(ds2) is the stored object:",
    identical(loo(ds2), ds2$criteria$loo), "\n")
cat("DONE\n")
