# What the standard error of marginaleffects::predictions() holds, against
# frm_linpred(se.fit = TRUE) and fitted(), for each kind of term get_coef() does
# not hold. Decides which fits the classed warning of get_vcov() reaches.
#   Rscript dev/predread-mfx-se.R <library>
# Data seed 7411.
args <- commandArgs(TRUE)
.libPaths(c(args[1L], "C:/Users/adf44/source/r/rellib-r16", .libPaths()))
suppressMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(7411)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(1:12, n, TRUE)))
ug <- rnorm(12, 0, 0.8)
d$cnt <- rpois(n, exp(0.5 + 0.4 * d$x + ug[d$g]))
d$yg <- 1 + 0.5 * d$x + ug[d$g] + rnorm(n, 0, 0.6)
d$lg <- exp(0.5 + 0.4 * d$x + ug[d$g] + rnorm(n, 0, 0.6))
d$sm <- sin(2 * d$x) + rnorm(n, 0, 0.4)
d$cgp <- rpois(n, exp(0.5 + sin(2 * d$x)))
nd <- data.frame(x = c(-1, 0, 1, 4), z = 0, g = factor(c(1, 2, 3, 3), 1:12))
show <- function(lab, fit, nd, type = "response") {
  calls <- 0L
  h <- function(w) invokeRestart("muffleWarning")
  p <- withCallingHandlers(predictions(fit, newdata = nd, type = type),
                           warning = h)
  fl <- suppressWarnings(frm_linpred(fit, newdata = nd, type = type,
                                     se.fit = TRUE))
  ft <- suppressWarnings(fitted(fit, newdata = nd))
  cat(sprintf("%-24s pred se %s\n%-24s linpred se %s\n%-24s fitted se %s\n",
              lab, paste(sprintf("%.4f", p$std.error), collapse = " "), "",
              paste(sprintf("%.4f", fl$se.fit), collapse = " "), "",
              paste(sprintf("%.4f", ft[, "Est.Error"]), collapse = " ")))
}
fits <- list(
  poisson_g = suppressWarnings(frm(cnt ~ x + (1 | g), family = poisson(), data = d)),
  gaussian_g = suppressWarnings(frm(yg ~ x + (1 | g), data = d)),
  lognormal_g = suppressWarnings(frm(lg ~ x + (1 | g), family = lognormal(), data = d)),
  s = suppressWarnings(frm(sm ~ s(x), data = d)),
  gp = suppressWarnings(frm(sm ~ gp(x, k = 15), data = d)),
  poisson_plain = suppressWarnings(frm(cnt ~ x, family = poisson(), data = d))
)
for (k in names(fits)) show(k, fits[[k]], nd)
cat("\nre_formula = NA:\n")
for (k in c("poisson_g", "gaussian_g", "lognormal_g")) {
  p <- suppressWarnings(predictions(fits[[k]], newdata = nd, re_formula = NA))
  ft <- suppressWarnings(fitted(fits[[k]], newdata = nd, re_formula = NA))
  cat(sprintf("%-14s pred %s | fitted %s\n", k,
              paste(sprintf("%.4f", p$std.error), collapse = " "),
              paste(sprintf("%.4f", ft[, "Est.Error"]), collapse = " ")))
}
cat("\nslope of x, poisson with (1 | g), and the same with the group ignored:\n")
print(suppressWarnings(avg_slopes(fits$poisson_g, variables = "x"))[, c("estimate", "std.error")])
cat("get_vcov dots on a predictions(re_formula = NA) call:\n")
tr <- list()
trace_ns <- asNamespace("frmtmb")
orig <- get("get_vcov.frmtmb_fit", trace_ns)
cat(" formals:", names(formals(orig)), "\n")
wrapped <- function(model, vcov = NULL, ...) {
  tr[[length(tr) + 1L]] <<- names(list(...))
  orig(model, vcov, ...)
}
registerS3method("get_vcov", "frmtmb_fit", wrapped,
                 envir = asNamespace("marginaleffects"))
invisible(suppressWarnings(predictions(fits$poisson_g, newdata = nd, re_formula = NA)))
cat(" calls:", length(tr), " dots:", paste(unlist(tr), collapse = ","), "\n")
tr <- list()
invisible(suppressWarnings(avg_slopes(fits$poisson_g, variables = "x")))
cat(" avg_slopes calls:", length(tr), " dots:", paste(unlist(tr), collapse = ","), "\n")
