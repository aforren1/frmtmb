# Lane fixes, item 1: brms 2.23.0's emmeans() on transformed predictors
# at frmtmb's estimates (Fixed_param draws, brms_fixed_fit()), beside
# frmtmb's and glm()'s.
#   Rscript dev/fixes-emm-brms.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
suppressMessages({library(testthat); library(frmtmb); library(emmeans)})
cat("LIB", find.package("frmtmb"), "\n")
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = globalenv())
  brms_rename <- frmtmb:::brms_rename
}
set.seed(20260930)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n, 1, 2), time = runif(n, 1, 5),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$yc <- rpois(n, d$time * exp(0.3 + 0.4 * d$x + 0.2 * d$z))
fmt <- function(v) paste(format(v, digits = 9), collapse = " ")
run <- function(lab, fo, specs, ...) {
  cat("==", lab, "\n")
  fit <- frm(bf(fo), data = d, family = poisson())
  ref <- glm(fo, data = d, family = poisson())
  bb <- brms_fixed_fit(brms::bf(fo), poisson(), d, fit, ndraws = 4)
  out <- list(
    frmtmb = tryCatch(summary(emmeans(fit, specs, ...))$emmean,
                      error = function(e) conditionMessage(e)),
    brms = tryCatch(summary(emmeans(bb, specs, ...))$emmean,
                    error = function(e) conditionMessage(e)),
    glm = tryCatch(summary(emmeans(ref, specs, ...))$emmean,
                   error = function(e) conditionMessage(e)))
  for (k in names(out)) cat(sprintf("  %-7s %s\n", k, fmt(out[[k]])))
  g <- tryCatch(ref_grid(bb, ...)@grid, error = function(e) NULL)
  if (!is.null(g)) cat("  brms grid z:", fmt(unique(g$z)), "\n")
  invisible(NULL)
}
run("poly(z, 2) + f", yc ~ poly(z, 2) + f, ~ f)
run("poly(z, 2) + f, at z = c(-1, 2)", yc ~ poly(z, 2) + f, ~ z + f,
    at = list(z = c(-1, 2)))
run("log(abs(z) + 1) + f", yc ~ log(abs(z) + 1) + f, ~ f)
run("scale(z) + f", yc ~ scale(z) + f, ~ f)
run("scale(z) + f, at z = 3", yc ~ scale(z) + f, ~ f, at = list(z = 3))
run("scale(z) + f, at z = c(0, 4)", yc ~ scale(z) + f, ~ z + f,
    at = list(z = c(0, 4)))
run("scale(z) + f + offset", yc ~ scale(z) + f + offset(log(time)), ~ f)
run("poly(z, 2) + f, epred", yc ~ poly(z, 2) + f, ~ f, epred = TRUE)
