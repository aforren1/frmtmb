# Reviewer, punch round 1: B1, update() and refits of `.` fits.
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "after"
base <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "after") c("C:/Users/adf44/source/r/wt-formula2-lib",
                                base) else base)
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
fx <- function(f) {
  if (inherits(f, "ERR")) paste("ERR:", substr(f, 1, 150)) else
    paste(rownames(fixef(f)), collapse = " ")
}
set.seed(6)
n <- 120
d <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n),
                g = factor(rep(1:6, each = 20)))
d$y <- d$y + d$x1
fit <- frm(y ~ ., data = d)
cat("== the stored call, formula(), print()\n")
cat("fit$call:\n"); print(fit$call)
cat("class of call$formula:", class(fit$call$formula), "\n")
cat("formula(fit):", deparse1(formula(fit)$formula %||% formula(fit)), "\n")
print(fit)
cat("\n== update()\n")
d2 <- d; d2$extra <- rnorm(n)
cat("newdata with extra column:", fx(q(update(fit, newdata = d2))), "\n")
cat("logLik identical to refit on d:",
    identical(logLik(q(update(fit, newdata = d2))), logLik(fit)), "\n")
cat(". ~ . + I(x1^2):", fx(q(update(fit, . ~ . + I(x1^2)))), "\n")
cat("formula. = ~ . - x1:", fx(q(update(fit, formula. = ~ . - x1))), "\n")
uf <- q(update(fit, family = student()))
cat("family = student():", fx(uf), "|",
    if (!inherits(uf, "ERR")) family(uf)$family, "\n")
dd <- d; dd$yo <- cut(dd$y, c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
dd$y <- NULL
fo <- frm(yo ~ ., data = dd, family = cumulative())
ua <- q(update(fo, bf(~ ., family = acat())))
cat("bf(~ ., family = acat()):", fx(ua), "|",
    if (!inherits(ua, "ERR")) family(ua)$family, "\n")
dw <- d; dw$w <- rnorm(n)
fw <- frm(y ~ . - w, data = dw)
cat("y ~ . - w, fixef:", fx(fw), "\n")
cat("newdata lacking the unused w:", fx(q(update(fw, newdata = d))), "\n")
cat("newdata lacking x2 (used):", fx(q(update(fit, newdata = d[, -3]))),
    "\n")
cat("\n== refit from the stored call after a column is added to d\n")
dref <- d
fr <- frm(y ~ ., data = dref)
dref$late <- rnorm(n)
cat("update(fr):", fx(q(update(fr))), "\n")
cat("eval(fr$call):", fx(q(eval(fr$call))), "\n")
cat("drop1(fr):", paste(rownames(q(drop1(fr))), collapse = " "), "\n")
cat("\n== a fit without a dot keeps its call\n")
f0 <- frm(y ~ x1 + g, data = d)
print(f0$call)
cat("update(f0, newdata = d2):", fx(q(update(f0, newdata = d2))), "\n")
cat("\n== multivariate and nonlinear fits with a dot\n")
dm <- d[, c("y", "x1", "x2")]; dm$y2 <- rnorm(n)
fm <- q(frm(bf(y ~ .) + bf(y2 ~ .) + set_rescor(FALSE), data = dm))
cat("mv fixef:", fx(fm), "\n")
dm2 <- dm; dm2$extra <- rnorm(n)
um <- q(update(fm, newdata = dm2))
cat("mv update(newdata + extra):", fx(um), "\n")
fn <- q(frm(bf(y ~ a + b * x1, a ~ ., b ~ 1, nl = TRUE), data = dm))
cat("nl fixef:", fx(fn), "\n")
un <- q(update(fn, newdata = dm2))
cat("nl update(newdata + extra):", fx(un), "\n")
fp <- q(frm(bf(y ~ ., sigma ~ x1) + student(), data = d))
up <- q(update(fp, newdata = d2))
cat("bf(y ~ ., sigma ~ x1) + student() update:", fx(up), "|",
    if (!inherits(up, "ERR")) family(up)$family, "\n")
cat("family arg on a dot fit:",
    fx(q(update(frm(y ~ ., data = d, family = gaussian()), newdata = d2))),
    "\n")
cat("DONE\n")
