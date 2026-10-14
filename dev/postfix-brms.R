# Lane postfix: brms 2.23.0's behavior on every item of the lane.
#   Rscript dev/postfix-brms.R
# Output: dev/postfix-log/brms.txt (by redirection).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
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
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## (c) threshold = as a brm() argument, before any compile
set.seed(4)
od <- data.frame(x = rnorm(100))
od$y <- cut(od$x + rlogis(100), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
show("brm(threshold = 'equidistant', chains = 0)",
     brm(y ~ x, family = cumulative(), data = od,
         threshold = "equidistant", chains = 0))
show("brm(cumulative(threshold =), chains = 0)",
     brm(y ~ x, family = cumulative(threshold = "equidistant"), data = od,
         chains = 0))
sc <- stancode(y ~ x, family = cumulative(), data = od,
               threshold = "equidistant")
cat("stancode(threshold = ) has delta:", grepl("delta", sc), "\n")
sc2 <- stancode(y ~ x, family = cumulative(threshold = "equidistant"),
                data = od)
cat("stancode(cumulative(threshold = )) has delta:", grepl("delta", sc2),
    "\n")

## (a) brm_multiple() on the brms_missings vignette's model, and a
## mixed one; 2 imputations, 1 chain each
imp <- mice::mice(mice::nhanes, m = 2, print = FALSE, seed = 1)
fm <- quiet(brm_multiple(bmi ~ age * chl, data = imp, chains = 1,
                         iter = 600, warmup = 300, seed = 1, refresh = 0))
dl <- mice::complete(imp, "all")
for (j in seq_along(dl)) dl[[j]]$g <- factor(rep(1:5, length.out = 25))
fmm <- quiet(brm_multiple(bmi ~ age + chl + (1 | g), data = dl,
                          chains = 1, iter = 600, warmup = 300, seed = 1,
                          refresh = 0))
for (obj in c("fm", "fmm")) {
  x <- get(obj)
  cat("--", obj, " class:", class(x), "\n")
  show("fixef", fixef(x))
  show("coef", coef(x))
  show("ranef", ranef(x))
  show("VarCorr", VarCorr(x))
  show("vcov", vcov(x))
  show("fitted", fitted(x))
  show("fitted(newdata)", fitted(x, newdata = dl[[1]][1:3, ]))
  show("predict", predict(x))
  show("residuals", residuals(x))
  show("conditional_effects(method=predict)",
       conditional_effects(x, method = "posterior_predict"))
  show("emmeans", emmeans::emmeans(x, ~ age))
  show("logLik", logLik(x))
  show("nobs", nobs(x))
  show("pp_check", pp_check(x))
  show("loo", loo(x))
  show("waic", waic(x))
  show("formula", formula(x))
  show("family", family(x))
  show("summary(waic = TRUE)", summary(x, waic = TRUE))
}
# what brms combines: the pooled fixef is the mean over the combined
# draws, and the SD over them
print(fixef(fm))
print(fitted(fm, newdata = dl[[1]][1:3, ]))
print(predict(fm, newdata = dl[[1]][1:3, ], seed = 1))
print(VarCorr(fmm))
print(nobs(fm))

## (d) summary(waic = TRUE) on a single fit: ignored?
s1 <- capture.output(print(summary(fm)))
s2 <- capture.output(print(summary(fm, waic = TRUE)))
cat("summary(waic = TRUE) prints the same as summary():",
    identical(s1, s2), "\n")
cat("formals(summary.brmsfit):",
    paste(names(formals(brms:::summary.brmsfit)), collapse = ", "), "\n")

## (b) pp_check(resp = ) on a multivariate fit
set.seed(3)
dd <- data.frame(x = rnorm(60))
dd$y1 <- 1 + dd$x + rnorm(60)
dd$y2 <- -1 + 0.5 * dd$x + rnorm(60)
mv <- quiet(brm(bf(mvbind(y1, y2) ~ x) + set_rescor(FALSE), data = dd,
                chains = 1, iter = 600, warmup = 300, seed = 1,
                refresh = 0))
show("pp_check(mv, resp = 'y1')", pp_check(mv, resp = "y1"))
show("pp_check(mv)", pp_check(mv))
show("pp_check(mv, resp = 'y2', type = 'stat')",
     pp_check(mv, resp = "y2", type = "stat"))
show("pp_check(mv, resp = c('y1', 'y2'))",
     pp_check(mv, resp = c("y1", "y2")))
show("pp_check(mv, resp = 'nope')", pp_check(mv, resp = "nope"))
p <- pp_check(mv, resp = "y1", ndraws = 5)
cat("pp_check(mv, resp = 'y1') layers:",
    paste(vapply(p$layers, function(l) class(l$geom)[1], ""),
          collapse = ", "), "\n")
mvr <- quiet(brm(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE), data = dd,
                 chains = 1, iter = 600, warmup = 300, seed = 1,
                 refresh = 0))
show("pp_check(mv rescor, resp = 'y1')", pp_check(mvr, resp = "y1"))
show("pp_check(mv rescor)", pp_check(mvr))

## (e) conditional_effects(method = "posterior_predict") on a fit
set.seed(5)
pd <- data.frame(x = rnorm(80))
pd$y <- rpois(80, exp(0.5 + 0.4 * pd$x))
pf <- quiet(brm(y ~ x, family = poisson(), data = pd, chains = 1,
                iter = 600, warmup = 300, seed = 1, refresh = 0))
ce <- conditional_effects(pf, method = "posterior_predict")
print(head(ce$x[, c("x", "estimate__", "se__", "lower__", "upper__")]))
ce2 <- conditional_effects(pf, method = "predict")
cat("method = 'predict' is posterior_predict:",
    identical(ce, ce2) || isTRUE(all.equal(ce$x$x, ce2$x$x)), "\n")
show("summary(pf, waic = TRUE)", summary(pf, waic = TRUE))
