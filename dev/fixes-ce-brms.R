# Lane fixes, item 3: what brms 2.23.0's conditional_effects() holds and
# computes for trunc() and se() terms whose variables are not in
# `conditions`, at frmtmb's estimates (Fixed_param draws), beside
# frmtmb.
#   Rscript dev/fixes-ce-brms.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
suppressMessages({library(testthat); library(frmtmb)})
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
cat("rows:", nrow(d), " mean(lo) =", format(mean(d$lo), digits = 10),
    " mean(s) =", format(mean(d$s), digits = 10),
    " mean(y) =", format(mean(d$y), digits = 10),
    " min(y) =", format(min(d$y), digits = 10), "\n")
fmt <- function(v) paste(format(v, digits = 7), collapse = " ")
xs <- c(-1, 0, 1)
run <- function(lab, fo_frm, fo_brms, method) {
  cat("==", lab, "| method =", method, "\n")
  fit <- frm(fo_frm, data = d, family = gaussian())
  bb <- tryCatch(brms_fixed_fit(fo_brms, gaussian(), d, fit, ndraws = 4000),
                 error = function(e) e)
  if (inherits(bb, "error")) {
    cat("  brms fixed fit ERROR:", conditionMessage(bb), "\n")
    return(invisible())
  }
  pc <- brms:::prepare_conditions(bb, effects = "x")
  cat("  brms conditions:", paste(names(pc), format(unlist(pc[1, ]),
                                                    digits = 10),
                                  sep = "=", collapse = " "), "\n")
  ce_b <- brms::conditional_effects(bb, "x", method = method,
                                    int_conditions = list(x = xs))[[1]]
  cat("  brms  estimate:", fmt(ce_b$estimate__), "\n")
  cat("  brms  lower   :", fmt(ce_b$lower__), "\n")
  cat("  brms  upper   :", fmt(ce_b$upper__), "\n")
  m <- if (method == "posterior_epred") "posterior_epred" else
    "posterior_predict"
  ce_f <- tryCatch(suppressMessages(
    conditional_effects(fit, "x", method = m,
                        int_conditions = list(x = xs))[[1]]),
    error = function(e) e)
  if (inherits(ce_f, "error")) {
    cat("  frm ERROR:", conditionMessage(ce_f), "\n")
  } else {
    held <- intersect(c("lo", "s", "y"), names(ce_f))
    cat("  frm held:", paste(held, format(unlist(ce_f[1, held]),
                                          digits = 10), sep = "=",
                             collapse = " "), "\n")
    cat("  frm   estimate:", fmt(ce_f$estimate__), "\n")
    cat("  frm   lower   :", fmt(ce_f$lower__), "\n")
    cat("  frm   upper   :", fmt(ce_f$upper__), "\n")
  }
  invisible()
}
models <- list(
  list("trunc(lb = lo)", bf(y | trunc(lb = lo) ~ x),
       brms::bf(y | trunc(lb = lo) ~ x)),
  list("trunc(ub = 10)", bf(y | trunc(ub = 10) ~ x),
       brms::bf(y | trunc(ub = 10) ~ x)),
  list("trunc(lb = min(y) - 1)", bf(y | trunc(lb = min(y) - 1) ~ x),
       brms::bf(y | trunc(lb = min(y) - 1) ~ x)),
  list("se(s)", bf(y2 | se(s) ~ x), brms::bf(y2 | se(s) ~ x)),
  list("se(s, sigma = TRUE)", bf(y2 | se(s, sigma = TRUE) ~ x),
       brms::bf(y2 | se(s, sigma = TRUE) ~ x)))
for (m in models) {
  for (meth in c("posterior_epred", "posterior_predict")) {
    run(m[[1]], m[[2]], m[[3]], meth)
  }
}
