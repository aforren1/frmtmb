# Item 3: every core family's default link, the nl body y ~ b0 + b1 * x from
# frm()'s default starts, against the linear fit. Data of dev/equiv-nl-sweep.R
# (seed 20261029). Also the non-default links with a bound at 0.
lib <- Sys.getenv("NLFIX_LIB",
                  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(20261029)
n <- 300
x <- round(rnorm(n), 3)
pos <- exp(0.4 + 0.3 * x + rnorm(n, 0, 0.4))
d <- data.frame(x = x, y = pos)
d$cnt <- rpois(n, exp(0.3 + 0.5 * x)) + 1L
d$w <- runif(n)
d$zz <- x + 3
ll <- function(expr) {
  msgs <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    invokeRestart("muffleWarning")
  }, message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  v <- if (inherits(r, "error")) NA_real_ else as.numeric(logLik(r))
  list(v = v, err = if (inherits(r, "error")) conditionMessage(r) else "",
       msg = paste(msgs, collapse = " "))
}
rows <- list(
  list("Gamma()", "y", "b0 + b1 * x", "x"),
  list("Gamma(link = 'inverse')", "y", "b1 * zz + b0", "zz"),
  list("inverse.gaussian()", "y", "b0 + b1 * x", "x"),
  list("Gamma(link = 'identity')", "y", "b0 + b1 * x", "x"),
  list("Gamma(link = 'sqrt')", "y", "b0 + b1 * x", "x"),
  list("Gamma(link = 'log')", "y", "b0 + b1 * x", "x"),
  list("inverse.gaussian(link = 'log')", "y", "b0 + b1 * x", "x"),
  list("inverse.gaussian(link = 'inverse')", "y", "b0 + b1 * x", "x"),
  list("poisson(link = 'identity')", "cnt", "b0 + b1 * x", "x"),
  list("poisson(link = 'sqrt')", "cnt", "b0 + b1 * x", "x"),
  list("lognormal()", "y", "b0 + b1 * x", "x"),
  list("weibull()", "y", "b0 + b1 * x", "x"),
  list("exponential()", "y", "b0 + b1 * x", "x"),
  list("frechet()", "y", "b0 + b1 * x", "x"),
  list("hurdle_gamma()", "y", "b0 + b1 * x", "x"))
for (r in rows) {
  fam <- eval(parse(text = r[[1]]))
  fl <- stats::as.formula(paste(r[[2]], "~", r[[4]]))
  fb <- stats::as.formula(paste(r[[2]], "~", r[[3]]))
  a <- ll(frm(bf(fl), family = fam, data = d))
  b <- ll(frm(bf(fb, b0 + b1 ~ 1, nl = TRUE), family = fam, data = d))
  cat(sprintf("%-38s linear %10.4f nl %10.4f diff %9.2e %s\n", r[[1]], a$v, b$v,
              b$v - a$v,
              if (nzchar(b$err)) paste("ERR", substr(b$err, 1, 120)) else ""))
  if (nzchar(b$msg)) cat("    message:", substr(b$msg, 1, 170), "\n")
}
cat("\nundefined at 0: keeps the error asking for a start\n")
b <- ll(frm(bf(y ~ log(a) + b1 * x, a + b1 ~ 1, nl = TRUE), family = Gamma(),
            data = d))
cat(substr(b$err, 1, 300), "\n")
cat("\nb1 * x alone, x changes sign: no start is in the domain\n")
b <- ll(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = Gamma(), data = d))
cat(substr(b$err, 1, 200), "\n")
cat("\nthe user's start stands\n")
b <- ll(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = Gamma(),
            data = d, start = list(beta = c("b0_(Intercept)" = 0.4))))
cat(b$v, "msg:", b$msg, "\n")
cat("\na located prior stands\n")
b <- ll(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = Gamma(),
            data = d, prior = set_prior("normal(0.4, 1)", nlpar = "b0")))
cat(b$v, "msg:", substr(b$msg, 1, 140), "\n")
cat("\na multivariate model moves only the body that needs it\n")
d$g2 <- x + rnorm(n)
b <- ll(frm(mvbf(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE, family = Gamma()),
                 bf(g2 ~ c0 + c1 * x, c0 + c1 ~ 1, nl = TRUE,
                    family = gaussian())), data = d))
cat(b$v, substr(b$err, 1, 200), "msg:", b$msg, "\n")
a1 <- ll(frm(bf(y ~ x), family = Gamma(), data = d))
a2 <- ll(frm(bf(g2 ~ x), family = gaussian(), data = d))
cat("sum of the two linear fits:", a1$v + a2$v, "\n")
print(par_template(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), data = d,
                   family = Gamma())$beta)
