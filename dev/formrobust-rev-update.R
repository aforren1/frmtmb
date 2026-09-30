# Reviewer: update() with a complete formula (claim 4). Data seed 41.
# REVLIB="" for the base arm.
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  w <- character(0); m <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) {
    m <<- c(m, trimws(conditionMessage(c))); invokeRestart("muffleMessage")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) cat(sprintf("[%s] ERROR: %s\n", label, r))
  if (length(w)) cat(sprintf("[%s] WARN: %s\n", label, unique(w)))
  if (length(m)) cat(sprintf("[%s] MSG: %s\n", label, unique(m)))
  invisible(if (inherits(r, "err")) NULL else r)
}
show <- function(label, fit) {
  if (is.null(fit)) return(invisible())
  cat(sprintf("[%s] coefs: %s | family: %s\n", label,
              paste(rownames(fixef(fit)), collapse = " "),
              fit$spec$responses[[1]]$family$family))
}
set.seed(41)
n <- 150
d <- data.frame(x = rnorm(n), x2 = rnorm(n), z = rnorm(n), w = rnorm(n),
                Age = runif(n, 1, 5))
d$y <- rnorm(n, 0.5 * d$x, exp(0.2 * d$z))
d$y2 <- rnorm(n, -0.3 * d$x)
d$count <- rpois(n, exp(0.5 + 0.2 * d$Age))
d$yp <- rpois(n, exp(0.3 * d$x))

f0 <- frm(bf(y ~ x, sigma ~ z), data = d)
show("orig", f0)
show("y ~ x2", tr("u1", update(f0, y ~ x2)))
show("bf(y ~ x2)", tr("u2", update(f0, bf(y ~ x2))))
show("bf(y ~ x2, sigma ~ w)", tr("u3", update(f0, bf(y ~ x2, sigma ~ w))))
show("bf(y ~ x2, sigma ~ 1)", tr("u4", update(f0, bf(y ~ x2, sigma ~ 1))))
show("~ . + x2", tr("u5", update(f0, ~ . + x2)))
show("y ~ x2, family = student()", tr("u6", update(f0, y ~ x2,
                                                   family = student())))
show("bf(y ~ x2, family = student())",
     tr("u7", update(f0, bf(y ~ x2, family = student()))))
show("bf(y ~ x2, family = student()) + family arg",
     tr("u8", update(f0, bf(y ~ x2, family = student()),
                     family = gaussian())))
u <- tr("u9", update(f0, y ~ x2))
if (!is.null(u)) {
  cat("[call of y ~ x2] ", deparse1(u$call), "\n")
  cat("[formula of y ~ x2] "); print(formula(u))
  # refit consistency: logLik equals a direct fit of the pooled model
  dfit <- frm(bf(y ~ x2, sigma ~ z), data = d)
  cat("[y ~ x2 logLik == direct fit]", identical(logLik(u), logLik(dfit)),
      "\n")
  # a second update from the updated fit
  u2 <- tr("u9b", update(u, y ~ x))
  if (!is.null(u2)) cat("[update of update logLik == f0]",
                        identical(logLik(u2), logLik(f0)), "\n")
  u3 <- tr("u9c", update(u, newdata = d[1:100, ]))
  show("update(newdata) of pooled", u3)
  u4 <- tr("u9d", update(u, evaluate = FALSE))
  if (!is.null(u4)) cat("[evaluate = FALSE] class", class(u4), "\n")
}
cat("\n-- a dpar constant --\n")
fc <- frm(bf(y ~ x, sigma = 1), data = d)
show("sigma = 1 orig", fc)
show("sigma = 1, update y ~ x2", tr("c1", update(fc, y ~ x2)))
show("sigma = 1, update bf(y ~ x2, sigma ~ z)",
     tr("c2", update(fc, bf(y ~ x2, sigma ~ z))))

cat("\n-- nonlinear --\n")
fn <- frm(bf(count ~ exp(a) * Age^b, a ~ 1, b ~ 1, nl = TRUE), data = d,
          family = poisson(link = "identity"))
show("nl orig", fn)
show("nl plain formula", tr("n1", update(fn, count ~ exp(a) * Age^(b / 2))))
show("nl bf nl=TRUE, a ~ x", tr("n2", update(fn, bf(count ~ exp(a) * Age^b,
                                                   a ~ x, nl = TRUE))))
show("nl bf without nl", tr("n3", update(fn, bf(count ~ exp(a) * Age^b))))
show("nl to linear bf(count ~ Age)", tr("n4", update(fn, bf(count ~ Age))))
show("nl to linear count ~ Age", tr("n5", update(fn, count ~ Age)))
tr("bf(nl = TRUE) no pforms at frm()",
   frm(bf(count ~ a + b, nl = TRUE), data = d, family = poisson()))
tr("bf(nl = TRUE) no pforms", print(class(bf(count ~ a + b, nl = TRUE))))
show("bf(nl) + lf(a ~ 1) + lf(b ~ 1)",
     tr("n6", frm(bf(count ~ exp(a) * Age^b, nl = TRUE) + lf(a ~ 1) +
                    lf(b ~ 1), data = d,
                  family = poisson(link = "identity"))))

cat("\n-- multivariate --\n")
fm <- frm(bf(y ~ x, sigma ~ z) + bf(y2 ~ x), data = d)
show("mv orig", fm)
um <- tr("m1", update(fm, bf(y ~ x2) + bf(y2 ~ x)))
if (!is.null(um)) cat("[mv update] coefs:",
                      rownames(fixef(um)), "\n")
um <- tr("m2", update(fm, ~ . + x2))
if (!is.null(um)) cat("[mv ~ . + x2] coefs:", rownames(fixef(um)), "\n")
cat("\n-- univariate model updated to a multivariate one --\n")
um <- tr("m3", update(f0, bf(y ~ x) + bf(y2 ~ x)))
if (!is.null(um)) cat("[uni -> mv] coefs:", rownames(fixef(um)), "\n")
