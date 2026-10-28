## brms 2.23.0 on mo() with two values (Jmo = 1) and one value (Jmo = 0):
## stancode, standata, default_prior, a fit, and what its methods report.
## Also marginaleffects' re_formula on a brms fit (item 5).
##   Rscript dev/mointerop-brms-mo12.R > dev/mointerop-out/brms-mo12.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(brms)
  library(marginaleffects)
})
cat("brms", format(packageVersion("brms")), "marginaleffects",
    format(packageVersion("marginaleffects")), "\n")
show <- function(tag, expr) {
  cat("\n==", tag, "==\n")
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n")
  else print(r)
  invisible(r)
}
## seed 5, 200 rows, x in 1..2, z continuous, g in 8 groups
set.seed(5)
n <- 200
d <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n),
                g = factor(sample(1:8, n, TRUE)))
d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(8)[d$g] * 0.5 + rnorm(n)
d$xo <- factor(d$x, levels = 1:2, ordered = TRUE)
d$x1 <- 4
d$x01 <- d$x - 1L

sc <- stancode(y ~ mo(x) + z, data = d)
cat("\n== stancode, two values (mo-related lines) ==\n")
cat(grep("mo|simo|Jmo|Xmo|bsp|con_simo", strsplit(sc, "\n")[[1]],
         value = TRUE), sep = "\n")
sd2 <- standata(y ~ mo(x) + z, data = d)
cat("\n== standata two values: Jmo", sd2$Jmo, "Ksp", sd2$Ksp,
    "Xmo_1 table:\n")
print(table(sd2$Xmo_1))
cat("con_simo_1:", sd2$con_simo_1, "\n")
show("default_prior two values", default_prior(y ~ mo(x) + z, data = d))
show("default_prior two values, ordered factor",
     default_prior(y ~ mo(xo) + z, data = d))
sdo <- standata(y ~ mo(xo) + z, data = d)
cat("ordered: Jmo", sdo$Jmo, "Xmo_1", head(sdo$Xmo_1), "\n")

sc1 <- stancode(y ~ mo(x1) + z, data = d)
cat("\n== stancode, one value (mo-related lines) ==\n")
cat(grep("mo|simo|Jmo|Xmo|bsp", strsplit(sc1, "\n")[[1]], value = TRUE),
    sep = "\n")
sd1 <- standata(y ~ mo(x1) + z, data = d)
cat("\n== standata one value: Jmo", sd1$Jmo, "Xmo_1 unique",
    unique(sd1$Xmo_1), "con_simo_1 length", length(sd1$con_simo_1), "\n")
show("default_prior one value", default_prior(y ~ mo(x1) + z, data = d))
show("brm one value", {
  f1 <- brm(y ~ mo(x1) + z, data = d, chains = 1, iter = 400, seed = 1,
            refresh = 0)
  list(ndraws = posterior::ndraws(f1), variables = variables(f1))
})

## the two-value fit
fit <- brm(y ~ mo(x) + z, data = d, chains = 2, iter = 2000, seed = 3,
           refresh = 0)
show("summary fixed", summary(fit)$fixed)
show("variables", variables(fit))
dr <- posterior::as_draws_df(fit)
show("simo_mox1 draws summary",
     posterior::summarise_draws(posterior::subset_draws(
       posterior::as_draws(fit), variable = "simo_mox1")))
show("bsp_mox summary", posterior::summarise_draws(
  posterior::subset_draws(posterior::as_draws(fit), variable = "bsp_mox")))
show("prior_summary", prior_summary(fit))
ce <- conditional_effects(fit, effects = "x")
show("conditional_effects x", ce$x[, c("x", "z", "estimate__", "se__")])
ce2 <- conditional_effects(fit, effects = "z")
show("conditional_effects z: x held at", unique(ce2$z$x))
show("conditional_effects z:x moderator values",
     unique(conditional_effects(fit, effects = "z:x")[[1]]$x))
nd <- data.frame(x = c(1, 2), z = 0)
show("fitted newdata x = 1, 2, z = 0", fitted(fit, newdata = nd))
show("posterior_epred x = 3", posterior_epred(fit,
                                              newdata = data.frame(x = 3, z = 0)))
show("posterior_epred x = 0", posterior_epred(fit,
                                              newdata = data.frame(x = 0, z = 0)))
show("posterior_epred x = 1.5", posterior_epred(fit,
                                                newdata = data.frame(x = 1.5, z = 0)))
## the linear reading: fitted(2) - fitted(1) is bsp
e <- posterior_epred(fit, newdata = nd)
cat("mean of epred(x=2) - epred(x=1):", mean(e[, 2] - e[, 1]),
    " mean bsp_mox:", mean(dr$bsp_mox), "\n")
show("ordered factor fit, stancode equal to integer", identical(
  stancode(y ~ mo(xo) + z, data = d), sc))

## item 5: marginaleffects re_formula on a brms fit
fr <- brm(y ~ z + (1 | g), data = d, chains = 2, iter = 1000, seed = 3,
          refresh = 0)
ndg <- data.frame(z = c(-1, 0, 1), g = factor(c(1, 2, 3), levels = 1:8))
show("predictions re_formula default (NULL)",
     as.data.frame(predictions(fr, newdata = ndg))[, c("z", "g",
                                                         "estimate")])
show("predictions re_formula = NA", as.data.frame(
  predictions(fr, newdata = ndg, re_formula = NA))[, c("z", "g",
                                                       "estimate")])
show("predictions re_formula = ~ 0", as.data.frame(
  predictions(fr, newdata = ndg, re_formula = ~ 0))[, c("z", "g",
                                                        "estimate")])
show("predictions re_formula = ~ (1 | g)", as.data.frame(
  predictions(fr, newdata = ndg, re_formula = ~ (1 | g)))[, c("z", "g",
                                                              "estimate")])
show("brms fitted re_formula = NA", fitted(fr, newdata = ndg,
                                           re_formula = NA))
show("brms fitted re_formula = NULL", fitted(fr, newdata = ndg))
show("avg_comparisons re_formula = NA", avg_comparisons(
  fr, variables = "z", re_formula = NA))
show("avg_slopes re_formula = NA", avg_slopes(fr, variables = "z",
                                              re_formula = NA))
show("predictions re_formula = NA, no newdata (first 3)", head(as.data.frame(
  predictions(fr, re_formula = NA))[, c("z", "g", "estimate")], 3))
show("predictions allow_new_levels", as.data.frame(predictions(
  fr, newdata = data.frame(z = 0, g = "new"), allow_new_levels = TRUE))[,
    c("z", "estimate")])
show("predictions re_formula = 'x' (a string)", predictions(
  fr, newdata = ndg, re_formula = "x"))
show("marginaleffects supported args for brmsfit",
     tryCatch(marginaleffects:::sanitize_model_specific.brmsfit,
              error = function(e) "no sanitize_model_specific.brmsfit"))
