# Lane fixes, item 3: emmeans(epred = TRUE) on trunc() and se() fits,
# frmtmb beside brms at fixed parameters; and conditional_effects() on
# frmtmb.sample draws of a trunc(lb = lo) fit.
#   Rscript dev/fixes-ce-emm.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
suppressMessages({library(testthat); library(frmtmb); library(emmeans)})
cat("LIB", find.package("frmtmb"), "\n")
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = globalenv())
}
brms_rename <- frmtmb:::brms_rename
set.seed(20261005)
n <- 150
d <- data.frame(x = runif(n, -1, 1), lo = runif(n, -1.5, 0),
                s = runif(n, 0.3, 1))
d$y <- 1 + d$x + rnorm(n)
d <- d[d$y > d$lo & d$y < 10, ]
d$y2 <- 1 + d$x + rnorm(nrow(d), 0, d$s)
fmt <- function(v) paste(format(v, digits = 7), collapse = " ")
models <- list(
  list("trunc(lb = lo)", bf(y | trunc(lb = lo) ~ x),
       brms::bf(y | trunc(lb = lo) ~ x)),
  list("trunc(lb = min(y) - 1)", bf(y | trunc(lb = min(y) - 1) ~ x),
       brms::bf(y | trunc(lb = min(y) - 1) ~ x)),
  list("se(s)", bf(y2 | se(s) ~ x), brms::bf(y2 | se(s) ~ x)))
for (m in models) {
  cat("==", m[[1]], "\n")
  fit <- frm(m[[2]], data = d, family = gaussian())
  bb <- brms_fixed_fit(m[[3]], gaussian(), d, fit, ndraws = 4)
  for (ep in c(FALSE, TRUE)) {
    a <- tryCatch(summary(emmeans(fit, ~ x, epred = ep,
                                  at = list(x = c(-1, 1))))$emmean,
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    b <- tryCatch(summary(emmeans(bb, ~ x, epred = ep,
                                  at = list(x = c(-1, 1))))$emmean,
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("  epred = %s\n    frm : %s\n    brms: %s\n", ep, fmt(a),
                fmt(b)))
  }
  g <- tryCatch(ref_grid(bb, epred = TRUE)@grid, error = function(e) NULL)
  if (!is.null(g)) cat("  brms grid columns:", names(g), "|",
                       fmt(unlist(g[1, ])), "\n")
  g <- tryCatch(ref_grid(fit, epred = TRUE)@grid, error = function(e) NULL)
  if (!is.null(g)) cat("  frm  grid columns:", names(g), "|",
                       fmt(unlist(g[1, ])), "\n")
}
# conditional_effects() on draws
suppressMessages(library(frmtmb.sample))
cat("frmtmb.sample:", find.package("frmtmb.sample"), "\n")
fit <- frm(bf(y | trunc(lb = lo) ~ x), data = d, family = gaussian())
ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 1)
ce <- tryCatch(conditional_effects(ds, "x", int_conditions =
                                     list(x = c(-1, 0, 1)))[[1]],
               error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(ce)) cat("draws CE:", ce, "\n") else {
  cat("draws CE held lo:", unique(ce$lo), " estimate:",
      fmt(ce$estimate__), "\n")
}
