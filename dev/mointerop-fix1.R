## Fix round 1 of lane mointerop: the review's cases
## (dev/reviews/2026-10-27-mointerop.md), on the lane and on rellib-r14.
##   Rscript dev/mointerop-fix1.R <lane|base> > dev/mointerop-out/fix1-<arm>.txt
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mointerop-lib", base)
} else base)
suppressMessages({
  library(frmtmb)
  library(marginaleffects)
})
cat("lib:", dirname(find.package("frmtmb")), "\n")
show <- function(tag, expr) {
  cat("\n==", tag, "==\n")
  w <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n")
  else print(r)
  for (m in w) cat("WARNING:", substr(m, 1, 300), "\n")
  invisible(r)
}
q <- function(e) suppressWarnings(suppressMessages(e))
se_of <- function(fit, nd, ...) predictions(fit, newdata = nd, ...)$std.error

## 1. sandwiches on gap, one-value and boundary fits
set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
dg$g <- factor(rep(1:20, each = 10))
fg <- q(frm(y ~ mo(x), data = dg))
mg <- lm(y ~ factor(x), data = dg)
nd5 <- data.frame(x = 1:5)
ndo <- data.frame(x = c(1L, 2L, 5L))
for (v in list("HC0", "HC1", ~ g)) {
  lab <- if (inherits(v, "formula")) "~ g" else v
  show(paste("gap fit, vcov =", lab, "at x = 1..5"), se_of(fg, nd5, vcov = v))
  show(paste("lm(factor(x)), vcov =", lab, "at x = 1, 2, 5"),
       se_of(mg, ndo, vcov = v))
}
show("gap fit, vcov = vcov(full = TRUE)",
     se_of(fg, nd5, vcov = vcov(fg, full = TRUE)))
show("gap fit, vcov_cluster(~ g, full = TRUE) diag",
     diag(vcov_cluster(fg, ~ g, "CR1", full = TRUE)))

set.seed(5)
n <- 200
d1 <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n),
                 g = factor(sample(1:8, n, TRUE)))
d1$y <- 1 + 0.8 * d1$x + 0.5 * d1$z + rnorm(n)
d1$x1 <- 4L
f1 <- q(frm(y ~ mo(x1) + z, data = d1))
m1 <- lm(y ~ z, data = d1)
nd1 <- data.frame(x1 = 4L, z = c(0, 1))
for (v in list("HC0", "HC1", ~ g)) {
  lab <- if (inherits(v, "formula")) "~ g" else v
  show(paste("one-value fit, vcov =", lab), se_of(f1, nd1, vcov = v))
  show(paste("lm(y ~ z), vcov =", lab), se_of(m1, data.frame(z = c(0, 1)),
                                              vcov = v))
}

set.seed(31)
n <- 300
# the review's construction (dev/review-mointerop-vcov.R)
d <- data.frame(x = rnorm(n), g = factor(sample(1:30, n, TRUE)),
                w = runif(n, 0.5, 2), off = log(runif(n, 1, 3)))
d$ypos <- rgamma(n, shape = 2, rate = 2 / exp(0.5 + 0.3 * d$x))
d$ybin <- rbinom(n, 1, pnorm(0.2 + 0.6 * d$x))
d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + d$off))
d$y <- 1 + 0.5 * d$x + rnorm(n) * (1 + abs(d$x))
fr <- q(frm(y ~ x + (1 | g), data = d))
show("boundary RE fit, sd(Intercept) estimate", VarCorr(fr)$g$sd[, 1])
ndr <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(d$g)))
show("boundary RE fit, vcov = ~ g", se_of(fr, ndr, vcov = ~ g,
                                          re_formula = NA))
show("boundary RE fit, default", se_of(fr, ndr, re_formula = NA))

## 3. HC1 and stata against lm; REML full matrix
m0 <- lm(y ~ x, data = d)
f0 <- frm(y ~ x, data = d)
nd0 <- data.frame(x = c(-1, 0, 1))
for (v in c("HC1", "stata")) {
  show(paste("frm", v), se_of(f0, nd0, vcov = v))
  show(paste("lm", v), se_of(m0, nd0, vcov = v))
}
for (v in c("HC3", "robust", "HAC", "kenward-roger", "OPG")) {
  show(paste("frm", v), se_of(f0, nd0, vcov = v))
}
frl <- q(frm(y ~ x + (1 | g), data = d, REML = TRUE))
show("REML full matrix", get_vcov(frl, vcov = vcov(frl, full = TRUE)))

## 4. Gamma(log): observed against expected information
fgm <- frm(ypos ~ x, data = d, family = Gamma(link = "log"))
mgm <- glm(ypos ~ x, data = d, family = Gamma(link = "log"))
show("Gamma(log) HC0 slope SE, frmtmb / vcovHC", sqrt(
  get_vcov(fgm, vcov = "HC0")[2, 2] /
    sandwich::vcovHC(mgm, type = "HC0")[2, 2]))

## 2. one-value fit downstream
show("one-value logLik", logLik(f1))
show("lm(y ~ z) logLik", logLik(m1))
show("one-value AIC, lm AIC", c(AIC(f1), AIC(m1)))
f1z <- frm(y ~ z, data = d1)
show("anova(y ~ z, y ~ mo(x1) + z)", anova(f1z, f1))
set.seed(1)
show("one-value predict()", predict(f1, newdata = nd1, ndraws = 4000))
set.seed(1)
show("y ~ z predict()", predict(f1z, newdata = data.frame(z = c(0, 1)),
                                ndraws = 4000))
show("(mo(x1) | g) on one value", frm(y ~ mo(x1) + z + (mo(x1) | g),
                                      data = d1))

## 5. get_predict() arguments
set.seed(41)
n <- 300
dr <- data.frame(x = rnorm(n), g = factor(sample(1:10, n, TRUE)))
dr$y <- 1 + 0.5 * dr$x + rnorm(10)[dr$g] + rnorm(n)
fre <- frm(y ~ x + (1 | g), data = dr)
ndg <- data.frame(x = c(-1, 0, 1), g = factor(1:3, levels = 1:10))
for (a in c("re.form", "allow.new.levels", "sample_new_levels", "ndraws",
            "draw_ids", "zitype")) {
  args <- list(fre, newdata = ndg)
  args[[a]] <- switch(a, re.form = NA, allow.new.levels = TRUE,
                      sample_new_levels = "gaussian", ndraws = 10,
                      draw_ids = 1:3, zitype = "response")
  show(paste("predictions(", a, "=)"), q(do.call(predictions, args))$estimate)
}
fs <- frm(bf(y ~ x, sigma ~ x), data = dr)
show("dpar = sigma estimate", q(predictions(fs, newdata = nd0,
                                            dpar = "sigma"))$estimate)
show("fitted(dpar = sigma)", fitted(fs, newdata = nd0, dpar = "sigma")[, 1])

## 6. the degenerate (mo(x) | g) gap fit, seed 23
set.seed(23)
n <- 300
d3 <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                 g = factor(sample(1:10, n, TRUE)))
d3$y <- 1 + (0.4 + rnorm(10, 0, 0.2)[d3$g]) * d3$x + 0.5 * d3$z + rnorm(n)
f3 <- q(frm(y ~ mo(x) + z + (mo(x) | g), data = d3))
show("seed 23 get_vcov", get_vcov(f3))

## 7. the padding count
set.seed(51)
n <- 600
G <- 20
d5 <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                 g = factor(sample(seq_len(G), n, TRUE)))
d5$y <- 1 + rnorm(G, 0, 0.6)[d5$g] + 0.4 * (d5$x - 1) + 0.5 * d5$z +
  rnorm(n)
f5 <- q(frm(y ~ mo(x) + z + (1 | g), data = d5))
show("gap warning count, 4 rows of a 20-level factor", se_of(
  f5, data.frame(x = c(1, 2, 3, 5), z = 0, g = factor(1, levels = 1:G))))
