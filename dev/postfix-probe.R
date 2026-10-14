# Lane postfix: every item of the lane, called the way brms's users call
# it, on one build.
#   Rscript dev/postfix-probe.R <arm> [lib]
# arm "before" reads rellib-r8 alone; "after" puts <lib> in front.
# Output: one line per call, OK <class/dim> or ERROR <message>.
a <- commandArgs(TRUE)
arm <- a[1]
lib <- if (length(a) > 1L) a[2] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(getwd(), "dev/stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("ARM", arm, "frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
show <- function(label, expr) {
  r <- tryCatch(withCallingHandlers(
    {
      v <- expr
      d <- if (is.null(dim(v))) length(v) else paste(dim(v), collapse = "x")
      paste0("OK ", paste(class(v), collapse = "/"), " [", d, "]")
    },
    warning = function(w) {
      cat("  [warning] ", substr(conditionMessage(w), 1, 160), "\n",
          sep = "")
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")),
    error = function(e) paste("ERROR", gsub("\n", " ", conditionMessage(e))))
  cat(sprintf("%-46s %s\n", label, substr(r, 1, 400)))
}

## (a) frm_multiple(): the brms_missings vignette's model, and a mixed one
set.seed(1)
imp <- mice::mice(mice::nhanes, m = 5, print = FALSE, seed = 1)
fm <- suppressMessages(frm_multiple(bmi ~ age * chl, data = imp))
dl <- mice::complete(imp, "all")
set.seed(2)
for (j in seq_along(dl)) dl[[j]]$g <- factor(rep(1:5, length.out = 25))
fmm <- suppressWarnings(suppressMessages(
  frm_multiple(bmi ~ age + chl + (1 | g), data = dl)))
for (obj in c("fm", "fmm")) {
  x <- get(obj)
  cat("--", obj, "\n")
  show("fixef", fixef(x))
  show("summary", summary(x))
  show("coef", coef(x))
  show("ranef", ranef(x))
  show("VarCorr", VarCorr(x))
  show("vcov", vcov(x))
  show("fitted", fitted(x))
  show("predict", predict(x))
  show("residuals", residuals(x))
  show("conditional_effects", conditional_effects(x))
  show("conditional_effects(method=predict)",
       conditional_effects(x, method = "predict"))
  show("emmeans", emmeans::emmeans(x, ~ age))
  show("logLik", logLik(x))
  show("nobs", nobs(x))
  show("AIC", AIC(x))
  show("pp_check", pp_check(x))
  show("posterior_summary", posterior_summary(x))
  show("hypothesis", hypothesis(x, "chl = 0"))
  show("as_draws_df", posterior::as_draws_df(x))
  show("loo", loo(x))
  show("formula", formula(x))
  show("family", family(x))
  show("model.frame", model.frame(x))
  show("update", update(x, . ~ . - chl))
  show("plot", plot(x))
  show("posterior_epred", posterior_epred(x))
}

## (b) pp_check(resp = ) on a multivariate ML fit
set.seed(3)
dd <- data.frame(x = rnorm(60))
dd$y1 <- 1 + dd$x + rnorm(60)
dd$y2 <- -1 + 0.5 * dd$x + rnorm(60)
mv <- suppressMessages(frm(bf(mvbind(y1, y2) ~ x), data = dd))
show("pp_check(mv, resp = 'y1')", pp_check(mv, resp = "y1"))
show("pp_check(mv)", pp_check(mv))
show("pp_check(mv, resp = 'y2', type = 'stat')",
     pp_check(mv, resp = "y2", type = "stat"))
uv <- suppressMessages(frm(y1 ~ x, data = dd))
show("pp_check(uv)", pp_check(uv))
show("pp_check(uv, resp = 'y1')", pp_check(uv, resp = "y1"))

## (c) threshold = as an frm() argument
set.seed(4)
od <- data.frame(x = rnorm(100))
od$y <- cut(od$x + rlogis(100), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
show("frm(threshold = 'equidistant')",
     frm(y ~ x, family = cumulative(), data = od,
         threshold = "equidistant"))
show("frm(family = cumulative(threshold =))",
     frm(y ~ x, family = cumulative(threshold = "equidistant"), data = od))
show("frm(thres = 1)", frm(y ~ x, family = cumulative(), data = od,
                           thres = 1))

## (d) summary(waic = TRUE)
show("summary(uv, waic = TRUE)", summary(uv, waic = TRUE))
show("summary(uv, waic = FALSE)", summary(uv, waic = FALSE))
show("summary(uv)", summary(uv))

## (f) brmsfit-methods:314 and :955 on the ported fixtures
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
fit1 <- brms_fixture(1)
newdata <- data.frame(Age = 0, visit = paste0("a", 1:100), Trt = 0,
                      count = 20, patient = 1, Exp = 2, volume = 0)
show(":314 fitted(old_levels, ndraws = 10)",
     fitted(fit1, newdata = newdata, allow_new_levels = TRUE,
            sample_new_levels = "old_levels", ndraws = 10))
show(":317 fitted(gaussian, ndraws = 1)",
     fitted(fit1, newdata = newdata, allow_new_levels = TRUE,
            sample_new_levels = "gaussian", ndraws = 1))
show("fitted(old_levels) no ndraws",
     fitted(fit1, newdata = newdata, allow_new_levels = TRUE,
            sample_new_levels = "old_levels"))
fit2 <- brms_fixture(2)
show(":955 update(fit2, bf(count ~ a + b, nl))",
     update(fit2, formula. = bf(count ~ a + b, nl = TRUE)))

## (e) conditional_effects(method = "predict") on draws
set.seed(5)
pd <- data.frame(x = rnorm(80))
pd$y <- rpois(80, exp(0.5 + 0.4 * pd$x))
pf <- suppressMessages(frm(y ~ x, family = poisson(), data = pd))
dr <- suppressWarnings(suppressMessages(
  frm_sample(pf, chains = 1, iter = 600, warmup = 300, seed = 1,
             refresh = 0)))
show("ce(draws, method = 'predict')",
     conditional_effects(dr, method = "predict"))
show("ce(draws, method = 'posterior_predict')",
     conditional_effects(dr, method = "posterior_predict"))
show("ce(draws)", conditional_effects(dr))
show("ce(fit, method = 'predict')",
     conditional_effects(pf, method = "predict"))
mvd <- suppressWarnings(suppressMessages(
  frm_sample(mv, chains = 1, iter = 600, warmup = 300, seed = 1,
             refresh = 0)))
show("pp_check(mv draws, resp = 'y1')", pp_check(mvd, resp = "y1"))
show("summary(draws, waic = TRUE)", summary(dr, waic = TRUE))
