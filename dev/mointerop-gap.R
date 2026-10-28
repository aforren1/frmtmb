## marginaleffects standard errors on a mo() fit with an empty category,
## against the exact reference: lm() on factor(x) with the ML sigma,
## which a mo() term with an interior simplex reparameterizes.
##   Rscript dev/mointerop-gap.R <lane|base> > dev/mointerop-out/gap-<arm>.txt
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
    format(packageVersion("frmtmb")), "\n")
f7 <- function(x) paste(format(signif(as.numeric(x), 7)), collapse = " ")
run <- function(tag, expr) {
  w <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  cat(tag, if (inherits(r, "error")) paste("ERROR:", conditionMessage(r))
      else f7(r), "\n")
  for (m in w) cat("   WARNING:", substr(m, 1, 220), "\n")
  invisible(r)
}
q <- function(e) suppressWarnings(suppressMessages(e))

## case 1: the mocode review's gap data (seed 7, x in 1, 2, 5)
set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
fg <- q(frm(y ~ mo(x), data = dg))
m <- lm(y ~ factor(x), data = dg)
s <- sqrt(sum(resid(m)^2) / nrow(dg)) / summary(m)$sigma
ex <- predict(m, newdata = data.frame(x = c(1L, 2L, 5L)), se.fit = TRUE)
cat("case 1, seed 7, y ~ mo(x), x in 1, 2, 5\n")
cat("exact SE at x = 1, 2, 5:            ", f7(ex$se.fit * s), "\n")
run("fitted() Est.Error at x = 1..5:     ",
    fitted(fg, newdata = data.frame(x = 1:5))[, 2])
run("predictions() SE at x = 1..5:       ",
    predictions(fg, newdata = data.frame(x = 1:5))$std.error)
run("predictions() SE again (no warning):",
    predictions(fg, newdata = data.frame(x = 1:5))$std.error)
run("predictions() estimate at x = 1..5: ",
    predictions(fg, newdata = data.frame(x = 1:5))$estimate)
run("fitted() estimate at x = 1..5:      ",
    fitted(fg, newdata = data.frame(x = 1:5))[, 1])
run("get_coef names:", names(get_coef(fg)))
cat("get_coef names:", names(get_coef(fg)), "\n")
run("get_vcov finite:", all(is.finite(get_vcov(fg))))
V <- vcov(m) * s^2
run("exact 1 -> 5 contrast, SE:          ",
    c(coef(m)[3], sqrt(V[3, 3])))
r <- run("avg_comparisons(x = c(1, 5)):      ", {
  a <- avg_comparisons(fg, variables = list(x = c(1, 5)))
  c(a$estimate, a$std.error)
})
run("exact 1 -> 2 contrast, SE:          ", c(coef(m)[2], sqrt(V[2, 2])))
run("avg_comparisons(x = c(1, 2)):      ", {
  a <- avg_comparisons(fg, variables = list(x = c(1, 2)))
  c(a$estimate, a$std.error)
})
run("avg_comparisons(x = c(1, 3)) (gap):", {
  a <- avg_comparisons(fg, variables = list(x = c(1, 3)))
  c(a$estimate, a$std.error)
})
run("exact mean prediction over the data, SE:", {
  X <- model.matrix(m)
  a <- colMeans(X)
  c(sum(a * coef(m)), sqrt(drop(t(a) %*% V %*% a)))
})
run("avg_predictions() over the data:    ", {
  a <- avg_predictions(fg)
  c(a$estimate, a$std.error)
})

## case 2: a covariate and a gap at 4 (seed 11, n = 300, x in 1, 2, 3, 5, 6)
set.seed(11)
n <- 300
d2 <- data.frame(x = sample(c(1L, 2L, 3L, 5L, 6L), n, replace = TRUE),
                 z = rnorm(n))
d2$y <- 1 + c(0, 0.5, 1.4, 1.7, 2.6, 3.0)[d2$x] + 0.5 * d2$z + rnorm(n)
f2 <- q(frm(y ~ mo(x) + z, data = d2))
m2 <- lm(y ~ factor(x) + z, data = d2)
s2 <- sqrt(sum(resid(m2)^2) / n) / summary(m2)$sigma
nd <- data.frame(x = 1:6, z = 0)
nd_obs <- nd[nd$x != 4, ]
cat("\ncase 2, seed 11, y ~ mo(x) + z, x in 1, 2, 3, 5, 6\n")
cat("logLik frm", format(as.numeric(logLik(f2)), digits = 12), "lm",
    format(as.numeric(logLik(m2)), digits = 12), "\n")
run("exact SE at x = 1, 2, 3, 5, 6, z = 0:",
    predict(m2, newdata = nd_obs, se.fit = TRUE)$se.fit * s2)
run("fitted() Est.Error at x = 1..6:      ",
    fitted(f2, newdata = nd)[, 2])
run("predictions() SE at x = 1..6:        ",
    predictions(f2, newdata = nd)$std.error)
run("avg_slopes(z) SE (exact next line):  ", avg_slopes(
  f2, variables = "z")$std.error)
run("exact z SE:                          ",
    sqrt(diag(vcov(m2)))["z"] * s2)
V2 <- vcov(m2) * s2^2
run("exact 2 -> 6, SE:                    ",
    c(coef(m2)["factor(x)6"] - coef(m2)["factor(x)2"],
      sqrt(V2["factor(x)6", "factor(x)6"] + V2["factor(x)2", "factor(x)2"] -
             2 * V2["factor(x)6", "factor(x)2"])))
run("avg_comparisons(x = c(2, 6)):        ", {
  a <- avg_comparisons(f2, variables = list(x = c(2, 6)))
  c(a$estimate, a$std.error)
})
run("hypotheses() SE of get_coef rows:    ", hypotheses(f2)$std.error)
cat("hypotheses() terms:", hypotheses(f2)$term, "\n")
run("vcov = 2 * get_vcov ratio at x = 1..6:", {
  V <- get_vcov(f2)
  predictions(f2, newdata = nd, vcov = 2 * V)$std.error /
    predictions(f2, newdata = nd)$std.error
})
