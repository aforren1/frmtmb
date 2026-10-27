## Fix round 2 check: marginaleffects, influence and insight on mo()
## fits, against the exact factor(k) reference of dev/review-mocode-se.R
## (lm with the ML sigma). Seeds 5, 11 and 7 as there.
##   Rscript dev/review-mocode-me.R
.libPaths(c("C:/Users/adf44/source/r/wt-mocode-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb lib:", dirname(system.file(package = "frmtmb")),
    " marginaleffects", format(packageVersion("marginaleffects")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
f7 <- function(x) paste(format(signif(as.numeric(x), 7)), collapse = " ")
me <- function(lab, expr) {
  w <- character(0)
  r <- withCallingHandlers(tryCatch(suppressMessages(expr),
                                    error = function(e) e),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  cat(sprintf("%-34s %s\n", lab, if (inherits(r, "error"))
    paste("ERROR:", substr(conditionMessage(r), 1, 200)) else f7(r)))
  if (length(w)) cat(sprintf("%-34s   %d warning(s): %s\n", "", length(w),
                             substr(w[1L], 1, 110)))
}
ml <- function(m) m$df.residual / length(residuals(m))
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
fa <- q(frm(y ~ mo(k) + z, data = d))
m <- lm(y ~ factor(k) + z, data = d)
V <- vcov(m) * ml(m)
nd <- data.frame(k = 1:4, z = 0)
cat("== y ~ mo(k) + z\n")
me("exact predictions SE", predict(m, nd, se.fit = TRUE)$se.fit * sqrt(ml(m)))
me("predictions() SE", marginaleffects::predictions(fa, newdata = nd)$std.error)
me("exact k 1->2 / 1->3 SE", sqrt(c(V[2, 2], V[3, 3])))
me("avg_comparisons k 1->2 SE", marginaleffects::avg_comparisons(
  fa, variables = list(k = c(1, 2)))$std.error)
me("avg_comparisons k 1->3 SE", marginaleffects::avg_comparisons(
  fa, variables = list(k = c(1, 3)))$std.error)
me("avg_comparisons k 2->3 est, SE", {
  r <- marginaleffects::avg_comparisons(fa, variables = list(k = c(2, 3)))
  c(r$estimate, r$std.error)
})
me("exact k 2->3 est, SE", c(coef(m)[3] - coef(m)[2],
                             sqrt(V[3, 3] + V[2, 2] - 2 * V[2, 3])))
me("exact z slope est, SE", c(coef(m)["z"], sqrt(V["z", "z"])))
me("avg_slopes(z) est, SE", {
  r <- marginaleffects::avg_slopes(fa, variables = "z")
  c(r$estimate, r$std.error)
})
me("slopes(z) at k = 2 SE", marginaleffects::slopes(
  fa, variables = "z", newdata = data.frame(k = 2L, z = 0))$std.error)
me("hypotheses(b2 = 0) est, SE", {
  r <- marginaleffects::hypotheses(fa, hypothesis = "b2 = 0")
  c(r$estimate, r$std.error)
})
me("predictions hypothesis k3 - k2", {
  r <- marginaleffects::predictions(fa, newdata = data.frame(k = 2:3, z = 0),
                                    hypothesis = "b2 - b1 = 0")
  c(r$estimate, r$std.error)
})
cat("get_coef names:", names(marginaleffects::get_coef(fa)), "\n")
cf <- marginaleffects::get_coef(fa)
f2 <- marginaleffects::set_coef(fa, cf)
me("set_coef round trip identical", identical(fitted(f2, newdata = nd)[, 1],
                                              fitted(fa, newdata = nd)[, 1]))
cf2 <- cf
cf2[grep("^zeta", names(cf2))[1L]] <- cf2[grep("^zeta", names(cf2))[1L]] + 0.1
me("set_coef zeta moves k = 1..4",
   fitted(marginaleffects::set_coef(fa, cf2), newdata = nd)[, 1] -
     fitted(fa, newdata = nd)[, 1])
me("get_vcov dim, finite", {
  Vg <- marginaleffects::get_vcov(fa)
  c(dim(Vg), all(is.finite(Vg)), identical(rownames(Vg), names(cf)))
})

cat("\n== y ~ mo(k) * z, seed 11\n")
set.seed(11)
dc <- data.frame(k = sample(1:4, n, replace = TRUE))
dc$z <- rnorm(n)
dc$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[dc$k] +
  dc$z * c(0.2, 1, 1.5, 1.7)[dc$k] + rnorm(n)
fc <- q(frm(y ~ mo(k) * z, data = dc))
mc <- lm(y ~ factor(k) * z, data = dc)
ndc <- data.frame(k = 1:4, z = 1)
me("exact predictions SE",
   predict(mc, ndc, se.fit = TRUE)$se.fit * sqrt(ml(mc)))
me("predictions() SE", marginaleffects::predictions(fc, newdata = ndc)$std.error)
me("exact slope of z at k = 1..4 SE", {
  Vc <- vcov(mc) * ml(mc)
  L <- cbind(0, 0, 0, 0, 1, diag(1, 4)[, 2:4])
  sqrt(diag(L %*% Vc %*% t(L)))
})
me("slopes(z) at k = 1..4 SE", marginaleffects::slopes(
  fc, variables = "z", newdata = data.frame(k = 1:4, z = 0))$std.error)

cat("\n== y ~ mo(k) + (mo(k) | g), seed 5\n")
fre <- q(frm(y ~ mo(k) + (mo(k) | g), data = d))
me("fitted(re_formula = NA) SE", fitted(fre, newdata = data.frame(k = 1:4),
                                        re_formula = NA)[, 2])
me("predictions(re_formula = NA) SE", marginaleffects::predictions(
  fre, newdata = data.frame(k = 1:4), re_formula = NA)$std.error)
me("predictions(g = 1) SE", marginaleffects::predictions(
  fre, newdata = data.frame(k = 1:4, g = "1"))$std.error)
me("fitted(g = 1) SE", fitted(fre, newdata = data.frame(k = 1:4,
                                                         g = "1"))[, 2])

cat("\n== influence and insight, y ~ mo(k) + z on 60 rows\n")
d6 <- d[1:60, ]
f6 <- q(frm(y ~ mo(k) + z, data = d6))
inf <- q(influence(f6))
cat("dfbeta colnames:", colnames(dfbeta(inf)), "\n")
cat("dfbetas colnames:", colnames(dfbetas(inf)), "\n")
me("cooks.distance head", head(cooks.distance(inf), 3))
cat("insight get_parameters:",
    paste(insight::get_parameters(f6)[[1]], collapse = " "), "\n")
cat("insight get_varcov dim:", dim(insight::get_varcov(f6)), "\n")

cat("\n== gap fit, seed 7\n")
set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
fg <- q(frm(y ~ mo(x), data = dg))
me("predictions() SE x = 1..5", marginaleffects::predictions(
  fg, newdata = data.frame(x = 1:5))$std.error)
me("again (warning once?)", marginaleffects::predictions(
  fg, newdata = data.frame(x = 1:5))$std.error)
me("fitted() SE x = 1..5", fitted(fg, newdata = data.frame(x = 1:5))[, 2])
