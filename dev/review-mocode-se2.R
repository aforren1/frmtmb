## Fix round 1 check: does the NaN at a gap category reach every
## consumer that reports an SE; is the in-sample SE of (mo(k) | g) and
## mo(k) * f the newdata SE on the same rows; (mo(k) | g) at known and
## new levels. Seeds 5 and 7 as in dev/review-mocode-se.R.
##   Rscript dev/review-mocode-se2.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", arm, "lib:", dirname(system.file(package = "frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
show <- function(lab, expr) {
  w <- character(0)
  r <- withCallingHandlers(tryCatch(suppressMessages(expr),
                                    error = function(e) e),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  cat(sprintf("%-40s %s\n", lab, if (inherits(r, "error"))
    paste("ERROR:", substr(conditionMessage(r), 1, 140)) else
      paste(format(signif(as.numeric(r), 7)), collapse = " ")))
  if (length(w)) cat(sprintf("%-40s   warning: %s\n", "",
                             substr(unique(w)[1L], 1, 100)))
}
set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
fg <- q(frm(y ~ mo(x), data = dg))
nd <- data.frame(x = 1:5)
cat("== gap fit, x = 1..5\n")
show("fitted() Est.Error", fitted(fg, newdata = nd)[, 2])
show("fitted() Q2.5", fitted(fg, newdata = nd)[, 3])
show("frm_linpred(se.fit) se", frm_linpred(fg, newdata = nd,
                                           se.fit = TRUE)$se.fit)
show("predict() Est.Error", predict(fg, newdata = nd, seed = 1)[, 2])
show("predict() Q2.5", predict(fg, newdata = nd, seed = 1)[, 3])
show("ce(x) se__", conditional_effects(fg, effects = "x")[[1]]$se__)
show("ce(x) lower__", conditional_effects(fg, effects = "x")[[1]]$lower__)
show("emmeans(~ x, at = 1:5) SE",
     summary(emmeans::emmeans(fg, ~ x, at = list(x = 1:5)))$SE)
show("frm_lp_basis sqrt(diag(A V A'))", {
  lb <- frm_lp_basis(fg, newdata = nd, re_formula = NA)
  sqrt(diag(as.matrix(lb$A %*% lb$V %*% t(lb$A))))
})
show("frm_lp_basis nonest flag", {
  lb <- frm_lp_basis(fg, newdata = nd, re_formula = NA)
  if (is.null(lb$nonest)) NA else as.numeric(lb$nonest)
})
show("marginaleffects predictions SE", marginaleffects::predictions(
  fg, newdata = nd)$std.error)

set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$f <- factor(sample(c("a", "b"), n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] +
  rnorm(8, sd = 0.4)[d$g] + d$k * rnorm(8, sd = 0.2)[d$g] + rnorm(n)
cat("\n== (mo(k) | g), seed 5\n")
fre <- q(frm(y ~ mo(k) + (mo(k) | g), data = d))
show("in-sample Est.Error rows 1:6", fitted(fre)[1:6, 2])
show("newdata = same rows Est.Error", fitted(fre, newdata = d[1:6, ])[, 2])
show("g = 1, k = 1..4", fitted(fre, newdata = data.frame(
  k = 1:4, g = "1"))[, 2])
show("re_formula = NA, k = 1..4", fitted(fre, newdata = data.frame(
  k = 1:4), re_formula = NA)[, 2])
show("new level, k = 1..4", fitted(fre, newdata = data.frame(
  k = 1:4, g = "new"), allow_new_levels = TRUE)[, 2])
show("estimates g = 1", fitted(fre, newdata = data.frame(
  k = 1:4, g = "1"))[, 1])
cat("\n== mo(k) * f, seed 5\n")
ff <- q(frm(y ~ mo(k) * f, data = d))
show("in-sample Est.Error rows 1:6", fitted(ff)[1:6, 2])
show("newdata = same rows Est.Error", fitted(ff, newdata = d[1:6, ])[, 2])
