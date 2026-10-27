## Fix round 1 check: do marginaleffects' SEs (its own J V J' over
## get_coef()) include the mo() simplex, and does frm_lp_basis() flag
## the gap rows in se_nonest? Exact reference as in
## dev/review-mocode-se.R (lm on factor(k), ML sigma). Seeds 5 and 7.
##   Rscript dev/review-mocode-se3.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", arm, "lib:", dirname(system.file(package = "frmtmb")),
    " marginaleffects", format(packageVersion("marginaleffects")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
f7 <- function(x) paste(format(signif(as.numeric(x), 7)), collapse = " ")
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
fa <- q(frm(y ~ mo(k) + z, data = d))
m <- lm(y ~ factor(k) + z, data = d)
s2 <- m$df.residual / n
nd <- data.frame(k = 1:4, z = 0)
ex <- predict(m, newdata = nd, se.fit = TRUE)$se.fit * sqrt(s2)
cat("exact SE k = 1..4, z = 0:          ", f7(ex), "\n")
cat("fitted() Est.Error:                 ", f7(q(fitted(fa, newdata = nd))[, 2]),
    "\n")
r <- tryCatch(q(marginaleffects::predictions(fa, newdata = nd)$std.error),
              error = function(e) conditionMessage(e))
cat("marginaleffects predictions SE:     ", f7(r), "\n")
V <- vcov(m) * s2
ex41 <- sqrt(V["factor(k)4", "factor(k)4"])
r <- tryCatch(q(marginaleffects::avg_comparisons(
  fa, variables = list(k = c(1, 4)))), error = function(e) e)
cat("avg_comparisons(k 1 -> 4): exact est", f7(coef(m)["factor(k)4"]),
    "SE", f7(ex41), "; marginaleffects",
    if (inherits(r, "error")) conditionMessage(r) else
      paste(f7(r$estimate), "SE", f7(r$std.error)), "\n")

set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
fg <- q(frm(y ~ mo(x), data = dg))
lb <- q(frm_lp_basis(fg, newdata = data.frame(x = 1:5), re_formula = NA))
cat("gap: frm_lp_basis se_nonest:", format(lb$se_nonest), "\n")
cat("gap: frm_lp_basis names:", paste(names(lb), collapse = " "), "\n")
cv <- tryCatch(q(frm_curve(fg, "x", newdata = data.frame(x = 1:5))),
               error = function(e) e)
cat("gap: frm_curve:", if (inherits(cv, "error")) conditionMessage(cv) else
  paste(capture.output(print(utils::head(as.data.frame(cv), 5))),
        collapse = "\n"), "\n")
