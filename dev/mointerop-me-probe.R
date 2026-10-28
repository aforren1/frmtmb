## marginaleffects on frmtmb fits: vcov = spellings and re_formula.
## Same data as dev/mointerop-brms-mo12.R (seed 5, 200 rows, g in 8).
##   Rscript dev/mointerop-me-probe.R <lane|base>
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
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "marginaleffects",
    format(packageVersion("marginaleffects")), "\n")
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
set.seed(5)
n <- 200
d <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n),
                g = factor(sample(1:8, n, TRUE)))
d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(8)[d$g] * 0.5 + rnorm(n)
fr <- frm(y ~ z + (1 | g), data = d)
ndg <- data.frame(z = c(-1, 0, 1), g = factor(c(1, 2, 3), levels = 1:8))
est <- function(p) as.data.frame(p)[, intersect(c("z", "g", "estimate",
                                                   "std.error"),
                                                 names(p))]
show("predictions default", est(predictions(fr, newdata = ndg)))
show("predictions re_formula = NA", est(predictions(fr, newdata = ndg,
                                                    re_formula = NA)))
show("predictions re_formula = ~ 0", est(predictions(fr, newdata = ndg,
                                                     re_formula = ~ 0)))
show("predictions re_formula = ~ (1 | g)", est(predictions(
  fr, newdata = ndg, re_formula = ~ (1 | g))))
show("fitted re_formula = NA", fitted(fr, newdata = ndg, re_formula = NA))
show("fitted re_formula = NULL", fitted(fr, newdata = ndg))
show("avg_slopes re_formula = NA", avg_slopes(fr, variables = "z",
                                              re_formula = NA))

## vcov = spellings, against lm and glm
fl <- frm(y ~ z, data = d)
ml <- lm(y ~ z, data = d)
nd <- data.frame(z = c(-1, 0, 1))
V <- get_vcov(fl)
cat("\nget_vcov(fl):\n"); print(V)
s0 <- predictions(fl, newdata = nd)$std.error
s0l <- predictions(ml, newdata = nd)$std.error
show("frm vcov = 2 V ratio", predictions(fl, newdata = nd,
                                          vcov = 2 * V)$std.error / s0)
show("lm vcov = 2 V ratio", predictions(ml, newdata = nd,
                                         vcov = 2 * vcov(ml))$std.error / s0l)
show("frm vcov = FALSE", est(predictions(fl, newdata = nd, vcov = FALSE)))
show("lm vcov = FALSE", est(predictions(ml, newdata = nd, vcov = FALSE)))
show("frm vcov = TRUE ratio", predictions(fl, newdata = nd,
                                           vcov = TRUE)$std.error / s0)
show("frm vcov = HC3", est(predictions(fl, newdata = nd, vcov = "HC3")))
show("lm vcov = HC3 ratio to default", predictions(
  ml, newdata = nd, vcov = "HC3")$std.error / s0l)
show("frm vcov = function", est(predictions(
  fl, newdata = nd, vcov = function(m) 2 * get_vcov(m))))
show("frm vcov = ~ g", est(predictions(fl, newdata = nd, vcov = ~ g)))
show("frm vcov = unnamed matrix", est(predictions(fl, newdata = nd,
                                                  vcov = unname(2 * V))))
show("frm vcov = wrong-size matrix", est(predictions(
  fl, newdata = nd, vcov = diag(2))))
show("frm hypotheses vcov = 2 V", hypotheses(fl, vcov = 2 * V))
show("frm avg_slopes vcov = 2 V", avg_slopes(fl, variables = "z",
                                             vcov = 2 * V))
show("frm avg_slopes default", avg_slopes(fl, variables = "z"))
## the sandwich spellings against lm and glm
show("frm vcov = HC0", predictions(fl, newdata = nd,
                                   vcov = "HC0")$std.error)
show("lm vcov = HC0", predictions(ml, newdata = nd, vcov = "HC0")$std.error)
show("sandwich::vcovHC(lm, HC0) direct", sqrt(diag(
  sandwich::vcovHC(ml, type = "HC0"))))
show("frm vcov = ~ g", predictions(fl, newdata = nd, vcov = ~ g)$std.error)
show("lm vcov = ~ g", predictions(ml, newdata = nd, vcov = ~ g)$std.error)
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$z))
gp <- glm(cnt ~ z, data = d, family = poisson())
fp <- frm(cnt ~ z, data = d, family = poisson())
show("glm poisson get_vcov(vcov = ~ g) SEs", sqrt(diag(get_vcov(gp, vcov = ~ g))))
show("frm poisson get_vcov(vcov = ~ g) SEs", sqrt(diag(get_vcov(fp, vcov = ~ g))))
show("glm poisson get_vcov(vcov = HC0) SEs", sqrt(diag(get_vcov(gp, vcov = "HC0"))))
show("frm poisson get_vcov(vcov = HC0) SEs", sqrt(diag(get_vcov(fp, vcov = "HC0"))))
show("frm RE fit vcov = HC0", predictions(fr, newdata = ndg,
                                          vcov = "HC0")$std.error)
show("frm RE fit vcov = ~ g", predictions(fr, newdata = ndg,
                                          vcov = ~ g)$std.error)
show("frm vcov = full outer matrix", predictions(
  fl, newdata = nd, vcov = vcov_cluster(fl, ~ g, "CR1",
                                        full = TRUE))$std.error)
show("frm vcov = 'unconditional'", predictions(fl, newdata = nd,
                                               vcov = "unconditional"))
show("frm vcov = 'kenward-roger'", est(predictions(fl, newdata = nd,
                                                   vcov = "kenward-roger")))
