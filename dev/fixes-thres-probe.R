# Lane fixes, item 4a: probe ord_thres_linpred() on the fit, and what
# frm_linpred() gives for a fixed disc.
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(62)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
f1 <- frm(y ~ x, family = cumulative(), data = d)
print(names(f1$frame$linpreds))
print(tryCatch(head(frm_linpred(f1, dpar = "disc", type = "response")),
               error = function(e) conditionMessage(e)))
f2 <- tryCatch(frm(bf(y ~ x, disc = 2), family = cumulative(), data = d),
               error = function(e) conditionMessage(e))
if (!is.character(f2)) {
  print(names(f2$frame$linpreds))
  print(tryCatch(head(frm_linpred(f2, dpar = "disc", type = "response")),
                 error = function(e) conditionMessage(e)))
  print(head(frmtmb:::ord_thres_linpred(f2), 2))
} else print(f2)
print(head(frmtmb:::ord_thres_linpred(f1), 2))
fg <- frm(yg | thres(gr = g) ~ x, family = cumulative(), data = d)
print(head(frmtmb:::ord_thres_linpred(fg), 3))
print(frmtmb:::ord_thres_linpred(fg, newdata = d[1:3, ]))
