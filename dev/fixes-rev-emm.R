# Reviewer of lane fixes, claim 1: emmeans() on transformed predictors,
# cases the lane did not run, against glm()/lm() through emmeans.
#   Rscript dev/fixes-rev-emm.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(emmeans); library(splines)})
cat("LIB", find.package("frmtmb"), "\n")
set.seed(7101)
n <- 150
d <- data.frame(x = rnorm(n), z = runif(n, 0.5, 4), w = runif(n, 1, 3),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$yc <- rpois(n, d$w * exp(0.2 + 0.3 * log(d$z) + 0.2 * (d$f == "b")))
d$yg <- 1 + 0.5 * d$z - 0.1 * d$z^2 + as.numeric(d$f) +
  0.3 * log(d$z) * (d$f == "c") + rnorm(n, 0, exp(-0.3 + 0.2 * d$z))
show <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-46s %s\n", lab, paste(r, collapse = " ")))
}
cmp <- function(lab, fo, family, spec, ...) {
  show(lab, {
    fit <- suppressWarnings(frm(bf(fo), data = d, family = family))
    ref <- if (family$family == "gaussian") lm(fo, data = d) else
      glm(fo, data = d, family = family)
    a <- summary(emmeans(fit, spec, ...))
    b <- summary(emmeans(ref, spec, ...))
    sprintf("n=%d maxrel=%.3g | %s", nrow(a),
            max(abs(a$emmean - b$emmean) / abs(b$emmean)),
            paste(format(a$emmean, digits = 8), collapse = " "))
  })
}
P <- poisson(); G <- gaussian()
# orthogonal poly: coefficients must come from training data, also at a
# grid value outside the training range
cmp("poly(z,2) raw=F at z=c(0,6,10)", yc ~ poly(z, 2) + f, P,
    c("z", "f"), at = list(z = c(0, 6, 10)))
cmp("poly(z,3,raw=TRUE)", yc ~ poly(z, 3, raw = TRUE) + f, P, "f")
cmp("ns(z,3)", yc ~ ns(z, 3) + f, P, "f")
cmp("ns(z,3) at z=c(1,5)", yc ~ ns(z, 3) + f, P, c("z", "f"),
    at = list(z = c(1, 5)))
cmp("bs(z,df=4) at z=c(0.6,3.9)", yc ~ bs(z, df = 4) + f, P, c("z", "f"),
    at = list(z = c(0.6, 3.9)))
cmp("I(z^2)+z", yc ~ z + I(z^2) + f, P, "f")
cmp("I(z^2)+z at z=7", yc ~ z + I(z^2) + f, P, "f", at = list(z = 7))
cmp("log(z)*f by f", yg ~ log(z) * f, G, ~ z | f, at = list(z = c(1, 2)))
cmp("f:log(z) emtrends-free", yg ~ f + f:log(z), G, "f")
cmp("scale(z)*f, at z=c(0.5,8)", yc ~ scale(z) * f, P, c("z", "f"),
    at = list(z = c(0.5, 8)))
cmp("poly(z,2)+f+offset(log(w))", yc ~ poly(z, 2) + f + offset(log(w)), P,
    "f")
cmp("poly + offset at w=2", yc ~ poly(z, 2) + f + offset(log(w)), P, "f",
    at = list(w = 2))
cmp("ns + offset type=response", yc ~ ns(z, 3) + f + offset(log(w)), P,
    "f", type = "response")
cmp("scale(x) + poly(z,2) two transforms", yc ~ scale(x) + poly(z, 2) + f,
    P, "f", at = list(x = 1, z = 2))

# a transformed predictor in sigma: the emmean of sigma (link scale) at
# z = 3 against the hand value from the fit's own coefficients
show("sigma ~ scale(z), dpar = sigma, at z = 3", {
  fit <- frm(bf(yg ~ f, sigma ~ scale(z)), data = d, family = G)
  fe <- fixef(fit)[, "Estimate"]
  hand <- fe[["sigma_Intercept"]] + fe[["sigma_scalez"]] *
    (3 - mean(d$z)) / sd(d$z)
  a <- summary(emmeans(fit, ~ z, dpar = "sigma", at = list(z = 3)))
  sprintf("emm %.10f hand %.10f rel %.3g", a$emmean, hand,
          abs(a$emmean - hand) / abs(hand))
})
show("sigma ~ poly(z,2), dpar = sigma, at z = c(1,3)", {
  fit <- frm(bf(yg ~ f, sigma ~ poly(z, 2)), data = d, family = G)
  fe <- fixef(fit)[, "Estimate"]
  P2 <- predict(poly(d$z, 2), c(1, 3))
  hand <- fe[["sigma_Intercept"]] + P2 %*% fe[grep("^sigma_poly", names(fe))]
  a <- summary(emmeans(fit, ~ z, dpar = "sigma", at = list(z = c(1, 3))))
  sprintf("emm %s hand %s maxrel %.3g", paste(format(a$emmean, digits = 9),
                                             collapse = ","),
          paste(format(hand, digits = 9), collapse = ","),
          max(abs(a$emmean - hand) / abs(hand)))
})
# mu's emmeans with sigma ~ scale(z) on the grid route (epred) and the
# mu-only design route: the grid must still hold z at mean(z)
show("epred, mu ~ poly(z,2)+f, sigma ~ scale(z)", {
  fit <- frm(bf(yg ~ poly(z, 2) + f, sigma ~ scale(z)), data = d,
             family = G)
  ref <- lm(yg ~ poly(z, 2) + f, data = d)
  rg <- emmeans::ref_grid(fit, epred = TRUE)
  a <- summary(emmeans(fit, "f", epred = TRUE))
  b <- summary(emmeans(fit, "f"))
  sprintf("grid z=%.8f mean(z)=%.8f epred %s design %s",
          unique(rg@grid$z), mean(d$z),
          paste(format(a$emmean, digits = 8), collapse = ","),
          paste(format(b$emmean, digits = 8), collapse = ","))
})
# random intercept, grid route at re_formula = NULL
show("poly(z,2) + (1|g) re_formula=NULL", {
  dd <- d; dd$g <- factor(rep(1:10, length.out = n))
  fit <- suppressWarnings(frm(bf(yc ~ poly(z, 2) + f + (1 | g)), data = dd,
                              family = P))
  a <- summary(emmeans(fit, "f", re_formula = NULL, at = list(z = 9)))
  b <- summary(emmeans(fit, "f", at = list(z = 9)))
  sprintf("reNULL %s | reNA %s", paste(format(a$emmean, digits = 8),
                                       collapse = ","),
          paste(format(b$emmean, digits = 8), collapse = ","))
})
# a predictor found from the formula environment (R's rule, the
# 2026-09-30 decision): z is not a column of data
show("poly(zz,2) with zz from the environment", {
  zz <- d$z
  d2 <- d[, c("yc", "f")]
  fit <- frm(bf(yc ~ poly(zz, 2) + f), data = d2, family = P)
  ref <- glm(yc ~ poly(zz, 2) + f, data = d2, family = P)
  b <- tryCatch(paste(format(summary(emmeans(ref, "f"))$emmean, digits = 8),
                      collapse = ","), error = function(e)
                        paste("glm ERROR", conditionMessage(e)))
  a <- tryCatch(paste(format(summary(emmeans(fit, "f"))$emmean, digits = 8),
                      collapse = ","), error = function(e)
                        paste("frm ERROR", conditionMessage(e)))
  paste("frm:", a, "| glm:", b)
})
show("scale(zz) env, no transform in data", {
  zz <- d$z
  d2 <- d[, c("yc", "f")]
  fit <- frm(bf(yc ~ zz + f), data = d2, family = P)
  a <- tryCatch(paste(format(summary(emmeans(fit, "f"))$emmean, digits = 8),
                      collapse = ","), error = function(e)
                        paste("frm ERROR", conditionMessage(e)))
  a
})
# a nonlinear body reading a transform of z in a by-variable spec
show("nl a + b*log(z), a ~ f, b ~ 1", {
  fit <- frm(bf(yg ~ a + b * log(z), a ~ 0 + f, b ~ 1, nl = TRUE),
             data = d, family = G)
  fe <- fixef(fit)[, "Estimate"]
  a <- summary(emmeans(fit, ~ f, at = list(z = 2.5), epred = TRUE))
  hand <- fe[grep("^a_f", names(fe))] + fe[["b_Intercept"]] * log(2.5)
  sprintf("maxrel %.3g", max(abs(a$emmean - hand) / abs(hand)))
})
# missing values in z: na.omit drops rows; the raw frame must align
show("poly(z,2) with NA rows in z and in f", {
  dd <- d; dd$z[c(3, 17)] <- NA; dd$f[40] <- NA
  fit <- frm(bf(yc ~ poly(z, 2) + f), data = dd, family = P)
  ref <- glm(yc ~ poly(z, 2) + f, data = dd, family = P)
  a <- summary(emmeans(fit, "f")); b <- summary(emmeans(ref, "f"))
  sprintf("maxrel %.3g gridz frm %.8f glm %.8f",
          max(abs(a$emmean - b$emmean) / abs(b$emmean)),
          unique(emmeans::ref_grid(fit)@grid$z),
          unique(emmeans::ref_grid(ref)@grid$z))
})
