# Item 4235: which terms drop1() and add1() list on fits with special
# terms. Usage: Rscript dev/silent-drop1-probe.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(11)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                k = sample(1:4, n, TRUE),
                g = factor(sample(letters[1:8], n, TRUE)))
d$xo <- d$x + rnorm(n, 0, 0.3)
d$sdx <- 0.3
d$y <- 1 + 0.6 * d$k + 0.4 * d$z + sin(2 * d$x) + rnorm(8)[d$g] +
  rnorm(n, 0, 0.5)
show <- function(lab, expr) {
  cat("==", lab, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) cat(" ", r, "\n") else print(r)
}
f_mo <- frm(y ~ mo(k) + z, data = d)
show("terms labels mo", attr(terms(f_mo), "term.labels"))
show("drop1 mo", drop1(f_mo, test = "Chisq"))
f_s <- frm(y ~ s(x) + z, data = d)
show("terms labels s", attr(terms(f_s), "term.labels"))
show("drop1 s", drop1(f_s, test = "Chisq"))
f_me <- frm(y ~ me(xo, sdx) + z, data = d)
show("terms labels me", attr(terms(f_me), "term.labels"))
show("drop1 me", drop1(f_me, test = "Chisq"))
f_re <- frm(y ~ x + z + (1 | g), data = d)
show("drop1 re", drop1(f_re, test = "Chisq"))
f_z <- frm(y ~ z, data = d)
show("add1 mo", add1(f_z, scope = ~ . + mo(k), test = "Chisq"))
show("add1 s", add1(f_z, scope = ~ . + s(x), test = "Chisq"))
show("add1 x", add1(f_z, scope = ~ . + x, test = "Chisq"))
show("add1 me", add1(f_z, scope = ~ . + me(xo, sdx), test = "Chisq"))
cat("logLik y ~ z:", format(as.numeric(logLik(f_z)), digits = 10),
    " y ~ me(xo, sdx) + z:", format(as.numeric(logLik(f_me)), digits = 10),
    "\n")
d2 <- d
d2$x[1:10] <- NA
show("univariate mi() predictor",
     frm(y ~ mi(x) + z, data = d2))
# what a hand refit gives, for the mo() row
f_k0 <- frm(y ~ z, data = d)
cat("hand LRT mo(k):",
    format(2 * (as.numeric(logLik(f_mo)) - as.numeric(logLik(f_k0))),
           digits = 10), "\n")
