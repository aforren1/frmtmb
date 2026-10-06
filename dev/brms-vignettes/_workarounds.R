# Is each workaround of the hand translation still needed?
#
# The vignette scripts record, at 0.42.0, a SPELLING or a substitute
# wherever the brms line did not run. This script runs the BRMS line
# itself, for every workaround a later release may have made
# unnecessary, so "no longer needed" is measured and not read off NEWS.
# Each row is labeled with the workaround it tests: the row passing
# means the workaround can go.
#
# Where frmtmb's refusal claims to follow brms, the brms side is RUN
# too, with brms::stancode() or brms::standata() (no compile), so the
# claim is checked against brms 2.23.0 rather than recalled.
#
#   BV_LIB=<libs> BV_OUT=<dir> BV_DIR=<this dir> Rscript _workarounds.R
#
# Write BV_OUT to its own directory: these rows are not vignette calls
# and must not enter the scoreboard.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = dirname(normalizePath(sub("^--file=", "", grep("^--file=",
    commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
), "_harness.R"))
bv_load()
bv_init("workarounds")
set.seed(1234)

## ---- data, as the vignettes build it ---------------------------------

b <- c(2, 0.75)
x <- rnorm(100)
y <- rnorm(100, mean = b[1] * exp(b[2] * x))
dat1 <- data.frame(x, y)
inv_logit <- function(x) 1 / (1 + exp(-x))
ability <- rnorm(300)
p <- 0.33 + 0.67 * inv_logit(ability)
answer <- ifelse(runif(300, 0, 1) < p, 1, 0)
dat_ir <- data.frame(ability, answer)
loss <- brms::loss
inhaler <- brms::inhaler
kidney <- brms::kidney
dat_smooth <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
data("nhanes", package = "mice", envir = globalenv())
imp <- mice::mice(nhanes, m = 3, print = FALSE, seed = 1234)
cbpp <- lme4::cbpp
phylo <- tryCatch(
  ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex"),
  error = function(e) NULL)
data_simple <- tryCatch(utils::read.table(
  "https://paul-buerkner.github.io/data/data_simple.txt", header = TRUE),
  error = function(e) NULL)
A <- if (!is.null(phylo)) ape::vcv.phylo(phylo) else NULL
cat("phylogenetic data downloaded:", !is.null(A) && !is.null(data_simple),
    "\n")
nlform <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
             ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE)

## ---- start = in place of the nonlinear priors -------------------------
#
# 0.42.0: every nonlinear model needed start = because frmtmb began at
# zero. brms_nonlinear and brms_multilevel write prior(..., nlpar = ).

bv("model", "WA start: brms_nonlinear fit1, no prior, no start", {
  frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), data = dat1)
}, "SPELLING", "start = list(beta = c(1, 0))")

bv("model", "WA start: brms_nonlinear fit1 with the vignette's prior1", {
  prior1 <- prior(normal(1, 2), nlpar = "b1") +
    prior(normal(0, 2), nlpar = "b2")
  frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), data = dat1,
      prior = prior1)
}, "SPELLING", "start = list(beta = c(1, 0))")

bv("model", "WA start: brms_nonlinear fit_loss with the vignette's priors", {
  frm(bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
         ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1, nl = TRUE),
      data = loss, family = gaussian(),
      prior = c(prior(normal(5000, 1000), nlpar = "ult"),
                prior(normal(1, 2), nlpar = "omega"),
                prior(normal(45, 10), nlpar = "theta")))
}, "SPELLING", "start = list(beta = c(5000, 1, 45))")

bv("model", "WA start: brms_nonlinear fit_loss, no prior, no start", {
  frm(bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
         ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1, nl = TRUE),
      data = loss, family = gaussian())
}, "SPELLING", "start = list(beta = c(5000, 1, 45))")

fit_loss1 <- bv("model", "WA start: brms_multilevel fit_loss1 with nlprior", {
  nlprior <- c(prior(normal(5000, 1000), nlpar = "ult"),
               prior(normal(1, 2), nlpar = "omega"),
               prior(normal(45, 10), nlpar = "theta"))
  frm(formula = nlform, data = loss, family = gaussian(), prior = nlprior)
}, "SPELLING",
"start = list(beta = c(5000, 1, 45)), and nlpar = becoming dpar =")

bv("model", "WA start: brms_nonlinear fit_ir2 with its prior", {
  frm(bf(answer ~ 0.33 + 0.67 * inv_logit(eta), eta ~ ability, nl = TRUE),
      data = dat_ir, family = bernoulli("identity"),
      prior = prior(normal(0, 5), nlpar = "eta"))
}, "SPELLING", "prior dropped, start = list(beta = c(0, 0))")

bv("model", "WA start: brms_nonlinear fit_ir3 with its priors", {
  frm(bf(answer ~ guess + (1 - guess) * inv_logit(eta),
         eta ~ 0 + ability, guess ~ 1, nl = TRUE),
      data = dat_ir, family = bernoulli("identity"),
      prior = c(prior(normal(0, 5), nlpar = "eta"),
                prior(beta(1, 1), nlpar = "guess", lb = 0, ub = 1)))
}, "SPELLING",
"start = list(beta = c(0, 0.3)), and the beta(1, 1) density dropped")

bv("model", "WA start: brms_multilevel SAMPLE fit_loss1 with nlprior", {
  nlprior <- c(prior(normal(5000, 1000), nlpar = "ult"),
               prior(normal(1, 2), nlpar = "omega"),
               prior(normal(45, 10), nlpar = "theta"))
  frm_sample(nlform, data = loss, family = gaussian(), prior = nlprior,
             chains = 1, iter = 400, warmup = 200, seed = 1234,
             refresh = 0)
}, "SPELLING", "set_prior(class = 'Intercept', dpar = 'ult') and start =")

## ---- conditional_effects() on nonlinear fits --------------------------

if (!is.null(fit_loss1)) {
  conditions <- data.frame(AY = unique(loss$AY))
  rownames(conditions) <- unique(loss$AY)
  bv("post", "WA band: conditional_effects(fit_loss1), default band", {
    conditional_effects(fit_loss1)
  }, "SPELLING", "band = 'boot'")
  me <- bv("post", "WA band: the vignette's method = 'predict' call", {
    conditional_effects(fit_loss1, conditions = conditions,
                        re_formula = NULL, method = "predict")
  }, "SPELLING", "method = 'predict' becoming band = 'boot'")
  bv("post", "WA facet: plot(me_loss, ncol = 5, points = TRUE)", {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off())
    plot(me, ncol = 5, points = TRUE)
  }, "SPELLING", "par(mfrow = c(2, 5)) set by the caller")
}

## ---- ordinal thresholds ------------------------------------------------
#
# brms 2.23.0 moved threshold = into the family. Its brm() still takes
# the JSS paper's threshold = argument into ..., so the brms side is
# checked for what it does with it.

bv("model", "WA threshold: sratio(threshold = 'equidistant') + cs prior", {
  frm(rating ~ period + carry + cs(treat) + (1 | subject),
      data = inhaler, family = sratio(threshold = "equidistant"),
      prior = set_prior("normal(-1,2)", coef = "treat"))
}, "SPELLING", "threshold = and the cs prior both dropped")

bv("post", "brms 2.23.0: brm(threshold = 'equidistant') as a brm() argument", {
  sc <- brms::stancode(rating ~ period + carry + cs(treat) + (1 | subject),
                       data = inhaler, family = brms::sratio(),
                       threshold = "equidistant")
  sc2 <- brms::stancode(rating ~ period + carry + cs(treat) + (1 | subject),
                        data = inhaler,
                        family = brms::sratio(threshold = "equidistant"))
  cat("delta in the program, argument form:", grepl("delta", sc),
      "  family form:", grepl("delta", sc2), "\n")
  c(argument = grepl("delta", sc), family = grepl("delta", sc2))
}, NA_character_, "")

bv("model", "WA threshold: cumulative(threshold = 'equidistant')", {
  d <- data.frame(x = rnorm(120))
  d$o <- cut(d$x + rnorm(120), 4, labels = FALSE)
  frm(o ~ x, data = d, family = cumulative(threshold = "equidistant"))
}, "SPELLING", "threshold = has no spelling")

## ---- surface and smooth displays ---------------------------------------

bv("post", "WA surface: conditional_effects(surface = TRUE) on t2()", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  f <- frm(y ~ t2(x1, x2), data = dat_smooth)
  plot(conditional_effects(f, surface = TRUE), ask = FALSE)
}, "SPELLING", "effects = 'x1:x2' curves at three values")

bv("post", "WA smooths: conditional_smooths()", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  f <- frm(bf(y ~ s(x1) + s(x2)), data = dat_smooth)
  plot(conditional_smooths(f), ask = FALSE)
}, "SPELLING", "conditional_effects() as the stand-in")

## ---- priors in brms's spelling -----------------------------------------

if (!is.null(A) && !is.null(data_simple)) {
  bv("model",
     "WA prior: brms_phylogenetics model_simple with its prior block", {
    frm(phen ~ cofactor + (1 | gr(phylo, cov = A)), data = data_simple,
        family = gaussian(), data2 = list(A = A),
        prior = c(prior(normal(0, 10), "b"),
                  prior(normal(0, 50), "Intercept"),
                  prior(student_t(3, 0, 20), "sd"),
                  prior(student_t(3, 0, 20), "sigma")))
  }, "SPELLING", "set_prior() strings, and the sigma prior on dpar = 'sigma'")
}

bv("model", "WA prior: brms_overview fit1 with its three set_prior() rows", {
  frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
      data = kidney, family = lognormal(),
      prior = c(set_prior("normal(0,5)", class = "b"),
                set_prior("cauchy(0,2)", class = "sd"),
                set_prior("lkj(2)", class = "cor")))
}, NA_character_, "")

bv("model", "WA prior: dirichlet on a mo() simplex (brms_monotonic fit4)", {
  d <- data.frame(income = factor(sample(c("below_20", "20_to_40",
                                           "40_to_100", "greater_100"),
                                         100, TRUE),
                                  levels = c("below_20", "20_to_40",
                                             "40_to_100", "greater_100"),
                                  ordered = TRUE))
  d$ls <- as.numeric(d$income) + rnorm(100)
  frm(ls ~ mo(income), data = d,
      prior = prior(dirichlet(c(2, 1, 1)), class = "simo",
                    coef = "moincome1"))
}, "SPELLING", "no spelling; fit4 is fit1 again")

## ---- mixtures ------------------------------------------------------------

bv("model", "WA mixture: bf(y ~ 1, mu1 ~ x1, mu2 ~ x1)", {
  frm(bf(y ~ 1, mu1 ~ x1, mu2 ~ x1), data = dat_smooth,
      family = mixture(gaussian(), gaussian()))
}, "SPELLING", "y ~ x1 as the main formula, only mu2 named")

bv("post", "brms 2.23.0: bf(y ~ 1, mu1 ~ x1, mu2 ~ x1) in stancode()", {
  sc <- brms::stancode(brms::bf(y ~ 1, mu1 ~ x1, mu2 ~ x1),
                       data = dat_smooth,
                       family = brms::mixture(stats::gaussian(),
                                              stats::gaussian()))
  cat("brms builds the program:", is.character(sc), "\n")
  "built"
}, NA_character_, "")

fit_mix2 <- bv("model", "frmtmb spelling of the mixture, for pp_mixture()", {
  frm(bf(y ~ x1, mu2 ~ x1), data = dat_smooth,
      family = mixture(gaussian(), gaussian()))
}, NA_character_, "")

bv("post", "WA mixture: pp_mixture(fit) on the ML fit", {
  utils::head(pp_mixture(fit_mix2))
}, "SPELLING", "mixture_probs() as the point-fit substitute")

## ---- custom families -----------------------------------------------------

bv("post", "WA custom: brms's custom_family() call verbatim", {
  custom_family(
    "beta_binomial2", dpars = c("mu", "phi"),
    links = c("logit", "log"),
    lb = c(0, 0), ub = c(1, NA),
    type = "int", vars = "vint1[n]"
  )
}, "SPELLING", "R lpdf in place of lb/ub/type/vars and the Stan block")

bb2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = list(mu = "logit", phi = "log"), type = "discrete",
  # the density brms_customfamilies.R uses; base lbeta() is not taped
  lpdf = function(y, dpars, aterms) {
    RTMBdist::dbetabinom(y, aterms$vint1, dpars$mu * dpars$phi,
                         (1 - dpars$mu) * dpars$phi, log = TRUE)
  })
bv("model", "WA custom: family = beta_binomial2 (brms's argument form)", {
  frm(incidence | vint(size) ~ period + (1 | herd), data = cbpp,
      family = bb2)
}, "SPELLING", "bf(...) + fam in place of family =")

## ---- the multiple-imputation and multivariate accessors ------------------

fit_imp1 <- frm_multiple(bmi ~ age * chl, data = imp)
bv("post", "WA translation: per-submodel fixef() at the 0.61.0 shape", {
  est <- vapply(fit_imp1$fits, function(f) fixef(f)[, "Estimate"],
                numeric(4))
  print(round(t(est), 3))
}, "SPELLING", "the 0.42.0 line read unlist(fixef(f)) as five numbers")

bv("post", "WA translation: hypothesis(fit_imp1, 'age = 0')", {
  hypothesis(fit_imp1, "age = 0")
}, NA_character_, "")

bv("post", "brms 2.23.0: a hypothesis with no relation", {
  brms::hypothesis(data.frame(age = rnorm(50)), "age")
}, NA_character_, "")

data("BTdata", package = "MCMCglmm", envir = globalenv())
fmv <- frm(bf(mvbind(tarsus, back) ~ sex + hatchdate + (1 | p | fosternest)
              + (1 | q | dam)) + set_rescor(TRUE), data = BTdata)
bv("post", "WA mv: pp_check(fit1, resp = 'tarsus') on the ML fit", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  pp_check(fmv, resp = "tarsus")
}, "SPELLING", "predict(resp = ) against the observed column")

bv("post", "WA mv: add_criterion(fit1, 'loo')", {
  add_criterion(fmv, "loo")
}, "SPELLING", "AIC() printed rather than attached")

## ---- the dots refusal of 0.58.0 -------------------------------------------
#
# brms's plot.brmsfit() has N, variable and regex; summary.brmsfit()
# has no waic and swallows it in its dots.

fk <- frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal())
bv("post", "dots: plot(fit, N = 2, ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fk, N = 2, ask = FALSE)
}, "BEHAVIOR", "N = absorbed by ...")
bv("post", "dots: plot(fit, variable = '^b', regex = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(fk, variable = "^b", regex = TRUE)
}, "BEHAVIOR", "variable = and regex = absorbed by ...")
bv("post", "dots: summary(fit, waic = TRUE)", {
  summary(fk, waic = TRUE)
}, "BEHAVIOR", "waic = absorbed by ...")
bv("post", "brms 2.23.0: formals of plot.brmsfit and summary.brmsfit", {
  pf <- names(formals(getS3method("plot", "brmsfit",
                                  envir = asNamespace("brms"))))
  sf <- names(formals(getS3method("summary", "brmsfit",
                                  envir = asNamespace("brms"))))
  cat("plot.brmsfit:", paste(pf, collapse = ", "), "\n")
  cat("summary.brmsfit:", paste(sf, collapse = ", "), "\n")
  c(N = "N" %in% pf, variable = "variable" %in% pf, waic = "waic" %in% sf)
}, NA_character_, "")
bv("post", "stancode(fit) on a frmtmb_fit, brms namespace loaded", {
  stancode(fk)
}, "BEHAVIOR", "a refusal naming frm_sample()")

bv_done()
