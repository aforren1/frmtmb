# Lane surface: one line per defect of the lane's brief, on one build.
#
#   Rscript dev/surface-repros.R lane|base > dev/surface-out/repros-<arm>.txt
#
# brms is loaded (never attached) for items 1 and 7 and attached at the
# end for item 6, which needs it on the search path after frmtmb.
arm <- commandArgs(TRUE)[1]
source("dev/surface-env.R")
surface_env(arm)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "| frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "| frmtmb.sample",
    as.character(packageVersion("frmtmb.sample")), "\n")
show <- surface_show
grDevices::pdf(NULL)

## 1. stancode(), standata(), pp_mixture() on a fit
set.seed(1)
d <- data.frame(x = rnorm(40)); d$y <- d$x + rnorm(40)
f <- frm(y ~ x, data = d)
show("1 stancode(fit), brms not loaded", stancode(f))
show("1 standata(fit), brms not loaded", standata(f))
show("1 pp_mixture(non-mixture fit)", pp_mixture(f))
set.seed(4)
dm <- data.frame(y = c(rnorm(60, -2), rnorm(60, 3)))
fm <- frm(y ~ 1, family = mixture(gaussian(), gaussian()), data = dm)
pm <- show("1 pp_mixture(gaussian mixture fit)", pp_mixture(fm))
if (!inherits(pm, "surface_err")) {
  cat("   dim", dim(pm), "| dimnames[[2]]", dimnames(pm)[[2]],
      "| [[3]]", dimnames(pm)[[3]], "\n")
  mp <- mixture_probs(fm)
  cat("   max |Estimate - mixture_probs()|",
      format(max(abs(pm[, "Estimate", ] - mp))), "\n")
}
set.seed(2)
n <- 300
do <- data.frame(x = rnorm(n))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * do$x, -1.5 * do$x) + rlogis(n)
do$y <- cut(lat, c(-Inf, -1, 0, 1, Inf), labels = FALSE)
fo <- show("1 fit ordinal mixture(cumulative, cumulative)",
           frm(y ~ x, family = mixture(cumulative(), cumulative()),
               data = do))
if (!inherits(fo, "surface_err")) {
  po <- show("1 pp_mixture(ordinal mixture fit)", pp_mixture(fo))
  if (!inherits(po, "surface_err")) cat("   dim", dim(po), "\n")
}
invisible(loadNamespace("brms"))
show("1 stancode(fit), brms loaded", stancode(f))
show("1 standata(fit), brms loaded", standata(f))
show("1 brms::stancode(fit)", brms::stancode(f))
show("1 brms::pp_mixture(fit)", brms::pp_mixture(fm))

## 2. plot() of a fit and brms's N, variable, regex
kidney <- brms::kidney
fk <- frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal())
show("2 plot(fit, N = 2, ask = FALSE)", plot(fk, N = 2, ask = FALSE))
show("2 plot(fit, variable = '^b', regex = TRUE)",
     plot(fk, variable = "^b", regex = TRUE))
show("2 plot(fit) control", plot(fk, ask = FALSE))

## 3. update() and a class-wide lkj prior the new formula cannot use
pr <- c(set_prior("normal(0,5)", class = "b"),
        set_prior("cauchy(0,2)", class = "sd"),
        set_prior("lkj(2)", class = "cor"))
f1 <- frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal(), prior = pr)
u <- show("3 update(fit, drop the correlation)",
          update(f1, formula. = ~ . - (1 + age | patient) + (1 | patient)))
if (!inherits(u, "surface_err")) print(prior_summary(u))
show("3 frm(direct, lkj on a model with no correlation)",
     frm(time | cens(censored) ~ age * sex + disease + (1 | patient),
         data = kidney, family = lognormal(),
         prior = set_prior("lkj(2)", class = "cor")))

## 4. frm_multiple() pooled post-processing
data("nhanes", package = "mice")
imp <- mice::mice(nhanes, m = 3, print = FALSE, seed = 1)
fmu <- frm_multiple(bmi ~ age * chl, data = imp)
fx <- show("4 fixef(frm_multiple)", fixef(fmu))
if (!inherits(fx, "surface_err")) print(fx)
cat("   rownames(x$pooled):", rownames(fmu$pooled), "\n")
show("4 summary(frm_multiple)", capture.output(summary(fmu)))
show("4 plot(frm_multiple, variable = '^b', regex = TRUE)",
     plot(fmu, variable = "^b", regex = TRUE))
show("4 x$rhats", if (is.null(fmu$rhats)) stop("NULL") else fmu$rhats)
show("4 conditional_effects(frm_multiple, 'age:chl')",
     conditional_effects(fmu, "age:chl"))

## 7. default_prior() of an order = "none" mixture, sum-to-zero component
set.seed(3)
dp <- data.frame(x = rnorm(200))
dp$y <- sample(1:5, 200, TRUE)
fam7 <- mixture(cumulative(threshold = "sum_to_zero"), cumulative())
dpf <- show("7 default_prior(frmtmb, stz component)",
            default_prior(y ~ x, data = dp, family = fam7))
if (!inherits(dpf, "surface_err")) {
  dpf <- as.data.frame(dpf)
  cat("   frmtmb Intercept rows by dpar:",
      paste(names(table(dpf$dpar[dpf$class == "Intercept"])),
            table(dpf$dpar[dpf$class == "Intercept"]), sep = "=",
            collapse = " "), "\n")
}
dpb <- show("7 default_prior(brms, stz component)",
            brms::default_prior(y ~ x, data = dp,
                                family = brms::mixture(
                                  brms::cumulative(threshold = "sum_to_zero"),
                                  brms::cumulative())))
if (!inherits(dpb, "surface_err")) {
  dpb <- as.data.frame(dpb)
  cat("   brms Intercept rows by dpar:",
      paste(names(table(dpb$dpar[dpb$class == "Intercept"])),
            table(dpb$dpar[dpb$class == "Intercept"]), sep = "=",
            collapse = " "), "\n")
}
s1 <- show("7 default_prior(frmtmb, one stz family)",
           default_prior(y ~ x, data = dp,
                         family = cumulative(threshold = "sum_to_zero")))
if (!inherits(s1, "surface_err")) {
  cat("   frmtmb one-family Intercept rows:",
      sum(as.data.frame(s1)$class == "Intercept"), "\n")
}

## 8. fitted() on an ordinal fit whose disc has no fixed column
set.seed(11)
n <- 120
d8 <- data.frame(x = round(runif(n, 0, 5), 1), z = rnorm(n),
                 g = factor(sample(letters[1:8], n, TRUE)))
lat <- sin(d8$x) + 0.5 * d8$z + rlogis(n)
d8$y <- factor(cut(lat, quantile(lat, c(0, 0.3, 0.6, 1)),
                   include.lowest = TRUE, labels = FALSE), ordered = TRUE)
for (rhs in c("0 + gp(x, k = 6)", "0 + (1 | g)")) {
  fit8 <- suppressWarnings(frm(bf(as.formula(paste("y ~ z")),
                                  as.formula(paste("disc ~", rhs))),
                               data = d8, family = cumulative()))
  show(paste("8 fitted(), disc ~", rhs), fitted(fit8))
  show(paste("8 fitted(newdata), disc ~", rhs), fitted(fit8,
                                                       newdata = d8[1:5, ]))
  show(paste("8 predict(), disc ~", rhs), predict(fit8))
}

## 10. summary(waic =) and add_criterion()
show("10 summary(fit, waic = TRUE)", capture.output(summary(fk, waic = TRUE)))
show("10 add_criterion(fit, 'loo')", add_criterion(fk, "loo"))
show("10 add_criterion(fit, 'waic')", add_criterion(fk, "waic"))

## 6. (cs(1) | g) with brms attached after frmtmb; categorical message
set.seed(1)
d6 <- data.frame(x = rnorm(200), g = factor(sample(letters[1:6], 200, TRUE)))
d6$y <- sample(1:4, 200, TRUE)
show("6 (cs(1) | g), brms not attached",
     frm(y ~ x + (cs(1) | g), family = sratio(), data = d6))
d6$yc <- factor(sample(c("a", "b", "c"), 200, TRUE))
fc <- frm(yc ~ x, family = categorical(), data = d6)
show("6 categorical CE predict", conditional_effects(fc, "x",
                                                     method = "predict"))
suppressPackageStartupMessages(library(brms))
show("6 (cs(1) | g), brms attached after frmtmb",
     frm(y ~ x + (cs(1) | g), family = sratio(), data = d6))
show("6 (cs(1) | g) on cumulative, brms attached",
     frm(y ~ x + (cs(1) | g), family = cumulative(), data = d6))
show("6 brms: (cs(1) | g), sratio, stancode",
     brms::stancode(y ~ x + (cs(1) | g), family = brms::sratio(), data = d6))
cat("REPROS DONE\n")
