# The other readers of a me() predictor on new rows, frmtmb against brms at the
# same parameter values (data seed 41): emmeans, conditional_effects, frm_curve,
# marginaleffects' get_predict, residuals(newdata =), simulate(newdata =).
#   Rscript dev/predread-me-readers.R <lib> > dev/predread-log/me-readers-<tag>.txt
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r16"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r16",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = E)
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
tr <- function(lab, e) {
  cat("\n---", lab, "\n")
  tryCatch(print(quiet(e)), error = function(err) {
    cat("ERROR:", conditionMessage(err), "\n")
  })
  invisible(NULL)
}
set.seed(41)
n <- 120
d <- data.frame(z = rnorm(n), x = rnorm(n), sx = 0.4)
d$y <- 1 + 0.8 * d$x + 0.3 * d$z + rnorm(n, 0, 0.5)
nd <- data.frame(x = c(-1, 0, 1), sx = c(0.2, 0.4, 0.8), z = c(0, 1, -1))
fit <- quiet(frm(bf(y ~ me(x, sx) + z), data = d, family = gaussian()))
bf1 <- E$brms_fixed_fit(brms::bf(y ~ me(x, sx) + z), stats::gaussian(), d, fit,
                        ndraws = 4000)
tr("frmtmb emmeans ~ z", as.data.frame(emmeans::emmeans(
  fit, ~ z, at = list(z = c(0, 1)))))
tr("brms emmeans ~ z", as.data.frame(emmeans::emmeans(
  bf1, ~ z, at = list(z = c(0, 1)))))
tr("frmtmb CE x", utils::head(conditional_effects(fit, "x")[[1]][
  , c("x", "sx", "estimate__", "se__")], 3))
tr("frmtmb CE z", utils::head(conditional_effects(fit, "z")[[1]][
  , c("z", "x", "sx", "estimate__", "se__")], 3))
tr("brms CE z", utils::head(conditional_effects(bf1, "z")[[1]][
  , c("z", "x", "sx", "estimate__", "se__")], 3))
tr("frmtmb frm_curve", frm_curve(fit, newdata = nd, simultaneous = FALSE)[
  , c("x", ".estimate", ".se")])
tr("brms linpred at the same rows", {
  l <- brms::posterior_linpred(bf1, newdata = nd, re_formula = NA)
  rbind(mean = colMeans(l), sd = apply(l, 2, sd))
})
tr("frmtmb curve without sx", frm_curve(fit, newdata = nd[c("x", "z")],
                                        simultaneous = FALSE)$.estimate)
tr("frmtmb fitted without sx", fitted(fit, newdata = nd[c("x", "z")]))
tr("frmtmb sdx = 0", {
  nd0 <- nd
  nd0$sx <- 0
  fitted(fit, newdata = nd0)
})
tr("frmtmb get_predict link", marginaleffects::get_predict(
  fit, newdata = nd, type = "link"))
tr("frmtmb residuals newdata", {
  nd2 <- nd
  nd2$y <- c(0, 1, 2)
  residuals(fit, newdata = nd2)
})
tr("brms residuals newdata", {
  nd2 <- nd
  nd2$y <- c(0, 1, 2)
  residuals(bf1, newdata = nd2)
})
tr("frmtmb simulate newdata sd", {
  s <- simulate(fit, nsim = 4000, seed = 1, newdata = nd)
  sapply(s, identity) |> t() |> apply(2, sd) |> utils::head(3)
})
tr("frmtmb predict sd", apply(predict(fit, newdata = nd, summary = FALSE,
                                      ndraws = 4000,
                                      propagate_error = FALSE), 2, sd))
tr("brms predict sd", apply(brms::posterior_predict(bf1, newdata = nd), 2, sd))
tr("frmtmb hypothesis", hypothesis(fit, "mexsx > 0")$hypothesis$Estimate)
