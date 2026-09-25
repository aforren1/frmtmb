# Lane wt-predfix: reproduce each filed defect before touching it.
#   PREDFIX_ARM=base Rscript dev/predfix-repro.R > dev/predfix-log/repro-base.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
options(digits = 10)
say <- function(...) cat(..., "\n", sep = "")
try_msg <- function(expr) {
  withCallingHandlers(
    tryCatch({ v <- expr; list(ok = TRUE, value = v) },
             error = function(e) list(ok = FALSE,
                                      value = conditionMessage(e))),
    warning = function(w) {
      say("  warning: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    })
}

say("\n== item 1: fitted() on a multivariate fit")
set.seed(20260921)
n <- 150
d <- data.frame(x = rnorm(n), g = factor(rep(1:15, each = 10)))
e <- matrix(rnorm(2 * n), n, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2))
d$y1 <- 1 + 0.5 * d$x + e[, 1]
d$y2 <- -0.3 * d$x + e[, 2]
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x))
fr <- frm(bf(mvbind(y1, y2) ~ x) + gaussian() + set_rescor(TRUE), data = d)
r <- try_msg(fitted(fr))
say("  fitted(mv): ", if (r$ok) paste(dim(r$value), collapse = " x ") else
  paste("ERROR:", r$value))
r <- try_msg(fitted(fr, resp = "y2"))
say("  fitted(mv, resp = 'y2'): ", if (r$ok) paste(dim(r$value),
  collapse = " x ") else paste("ERROR:", r$value))
set.seed(1)
r <- try_msg(predict(fr, ndraws = 50))
say("  predict(mv): ", if (r$ok) paste(dim(r$value), collapse = " x ") else
  paste("ERROR:", r$value))

say("\n== item 2: the overflow fixture of dev/shapes-p2-measure.R")
fp <- frm(bf(cnt ~ x) + poisson(), data = d)
b <- fixef_by_dpar(fp)$mu
x_bad <- (709.5 - b[["(Intercept)"]]) / b[["x"]]
say("  x_bad = ", x_bad, " (data range ", paste(range(d$x), collapse = " "),
    ")")
set.seed(7)
r <- try_msg(predict(fp, newdata = data.frame(x = c(4, x_bad)), ndraws = 400))
print(r$value)

say("\n== item 4: a quadrature fit's scalar Est.Error")
set.seed(10)
db <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(60))
db$y <- rbinom(60, 1, plogis(db$x + rnorm(10)[db$g]))
fq <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = db,
                           quadrature = TRUE))
nd <- db[1:2, c("x", "g")]
r1 <- try_msg(fitted(fq, newdata = nd))
r2 <- try_msg(fitted(fq, newdata = nd, re_formula = NA))
say("  known level Est.Error: ",
    paste(format(r1$value[, "Est.Error"], digits = 6), collapse = " "))
say("  re_formula = NA Est.Error: ",
    paste(format(r2$value[, "Est.Error"], digits = 6), collapse = " "))
say("  finite-difference route (scale = 'linear'): ")
r3 <- try_msg(fitted(fq, newdata = nd, scale = "linear"))
print(r3$value)

say("\n== item 5: vcov() on a no-RE REML / profile fit")
set.seed(5)
d5 <- data.frame(x = rnorm(80))
d5$cnt <- rpois(80, exp(0.2 + 0.5 * d5$x))
d5$y <- 1 + d5$x + rnorm(80)
for (fam in c("gaussian", "poisson")) {
  for (mode in c("REML", "profile")) {
    fo <- if (fam == "gaussian") bf(y ~ x) + gaussian() else
      bf(cnt ~ x) + poisson()
    f <- if (mode == "REML") frm(fo, data = d5, REML = TRUE) else
      frm(fo, data = d5, control = frmtmb_control(profile = TRUE))
    r <- try_msg(vcov(f))
    s <- try_msg(summary(f))
    say("  ", fam, " ", mode, ": vcov ", if (r$ok) "answers" else
      paste("ERROR:", r$value), "; summary ", if (s$ok) "answers" else
      paste("ERROR:", s$value))
  }
}

say("\n== item 6: poisson y ~ 0 + x at covariate scale s, against glm")
set.seed(23)
n6 <- 250
xt <- rnorm(n6)
y6 <- rpois(n6, pmin(exp(0.5 + 0.4 * xt), 1e6))
for (s in c(1, 1e-2, 1e-4, 1e-6, 1e-8)) {
  dd <- data.frame(y = y6, xs = xt * s)
  gl <- as.numeric(logLik(glm(y ~ 0 + xs, family = poisson(), data = dd)))
  f0 <- withCallingHandlers(frm(bf(y ~ 0 + xs), family = poisson(), data = dd),
                            warning = function(w) {
                              say("  warning: ", conditionMessage(w))
                              invokeRestart("muffleWarning")
                            })
  say(sprintf("  scale %-6g glm %.6f default %.6f gap %.6f conv %d %s", s, gl,
              as.numeric(logLik(f0)), as.numeric(logLik(f0)) - gl,
              f0$opt$convergence, f0$opt$message))
}
