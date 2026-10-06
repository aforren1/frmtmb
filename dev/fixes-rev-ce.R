# Reviewer of lane fixes, claim 3: conditional_effects() on trunc() and
# se() without conditions, against brms 2.23.0 at frmtmb's estimates.
#   Rscript dev/fixes-rev-ce.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/fixes-rev-stan-cache"))
suppressMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = globalenv())
}
brms_rename <- frmtmb:::brms_rename
set.seed(4242)
n <- 200
d <- data.frame(x = runif(n, -1, 1), lo = runif(n, -2, 0.5),
                s = runif(n, 0.2, 1.2), g = factor(sample(letters[1:3], n,
                                                          TRUE)))
d$y <- 0.5 + d$x + rnorm(n)
d <- d[d$y > d$lo & d$y > d$x - 2, ]
d$y2 <- 0.5 + d$x + rnorm(nrow(d), 0, sqrt(0.5^2 + d$s^2))
cat("rows", nrow(d), "mean(lo)", format(mean(d$lo), digits = 10),
    "mean(s)", format(mean(d$s), digits = 10), "mean(y)",
    format(mean(d$y), digits = 10), "\n")
xs <- c(-1, 0, 1)
fmt <- function(v) paste(format(v, digits = 8), collapse = " ")
run <- function(lab, fo, bfo, cond = NULL) {
  cat("==", lab, if (!is.null(cond)) paste("conditions:",
                                           deparse1(cond)), "\n")
  fit <- frm(fo, data = d, family = gaussian())
  bb <- brms_fixed_fit(bfo, gaussian(), d, fit, ndraws = 4000)
  for (m in c("posterior_epred", "posterior_predict")) {
    set.seed(1)
    cb <- tryCatch(brms::conditional_effects(
      bb, "x", method = m, int_conditions = list(x = xs),
      conditions = cond)[[1]], error = function(e) conditionMessage(e))
    set.seed(1)
    # frmtmb refuses conditions = NULL (brms takes it), so it is left out
    fa <- list(fit, "x", method = m, int_conditions = list(x = xs),
               ndraws = 4000)
    if (!is.null(cond)) fa$conditions <- cond
    cf <- tryCatch(suppressMessages(do.call(conditional_effects, fa))[[1]],
                   error = function(e) conditionMessage(e))
    if (is.character(cb) || is.character(cf)) {
      cat("  ", m, "brms:", if (is.character(cb)) cb else "ok",
          "frm:", if (is.character(cf)) cf else "ok", "\n")
      next
    }
    held <- intersect(c("lo", "s", "y", "y2"), names(cf))
    cat("  ", m, "\n")
    cat("    held brms:", paste(held, format(unlist(cb[1, held]),
                                             digits = 10), collapse = " "),
        "\n    held frm :", paste(held, format(unlist(cf[1, held]),
                                               digits = 10), collapse = " "),
        "\n")
    cat("    est  brms:", fmt(cb$estimate__), "\n    est  frm :",
        fmt(cf$estimate__), "\n")
    cat("    low  brms:", fmt(cb$lower__), "\n    low  frm :",
        fmt(cf$lower__), "\n")
    cat("    upp  brms:", fmt(cb$upper__), "\n    upp  frm :",
        fmt(cf$upper__), "\n")
    if (m == "posterior_epred") {
      cat(sprintf("    est max rel diff %.3g\n",
                  max(abs(cb$estimate__ - cf$estimate__) /
                        abs(cb$estimate__))))
    }
  }
  # the predictive draws themselves must respect the bound: simulate at
  # the held grid and report the minimum against the bound
  invisible(fit)
}
run("trunc(lb = min(y) - 1)", bf(y | trunc(lb = min(y) - 1) ~ x),
    brms::bf(y | trunc(lb = min(y) - 1) ~ x))
run("trunc(lb = x - 2), bound on the effect",
    bf(y | trunc(lb = x - 2) ~ x), brms::bf(y | trunc(lb = x - 2) ~ x))
run("se(s, sigma = TRUE)", bf(y2 | se(s, sigma = TRUE) ~ x),
    brms::bf(y2 | se(s, sigma = TRUE) ~ x))
run("trunc(lb = lo), conditions lo = 0.25",
    bf(y | trunc(lb = lo) ~ x), brms::bf(y | trunc(lb = lo) ~ x),
    cond = data.frame(lo = 0.25))
run("se(s), conditions s = 2", bf(y2 | se(s) ~ x),
    brms::bf(y2 | se(s) ~ x), cond = data.frame(s = 2))
run("trunc(lb = lo) + g, factor held",
    bf(y | trunc(lb = lo) ~ x + g), brms::bf(y | trunc(lb = lo) ~ x + g))
run("weights(s) beside trunc(lb = lo)",
    bf(y | weights(s) + trunc(lb = lo) ~ x),
    brms::bf(y | weights(s) + trunc(lb = lo) ~ x))

# do the predictive draws respect the held bound? frmtmb's lower__ at a
# small probability, and the minimum of the draws
fit <- frm(bf(y | trunc(lb = lo) ~ x), data = d)
ce <- suppressMessages(conditional_effects(fit, "x", method = "predict",
                                           int_conditions = list(x = xs),
                                           ndraws = 20000, prob = 0.9999))[[1]]
cat("trunc(lb = lo) lower__ at prob 0.9999:", fmt(ce$lower__),
    " held bound", format(mean(d$lo), digits = 10), " all >= bound:",
    all(ce$lower__ >= mean(d$lo)), "\n")
fit <- frm(bf(y | trunc(lb = x - 2) ~ x), data = d)
ce <- suppressMessages(conditional_effects(fit, "x", method = "predict",
                                           int_conditions = list(x = xs),
                                           ndraws = 20000, prob = 0.9999))[[1]]
cat("trunc(lb = x - 2) lower__ at prob 0.9999:", fmt(ce$lower__),
    " bounds", fmt(xs - 2), " all >= bound:", all(ce$lower__ >= xs - 2), "\n")
# cens() on the same grid: no refusal before or after
d$cn <- ifelse(d$y > 2, "right", "none")
d$yc <- pmin(d$y, 2)
fit <- frm(bf(yc | cens(cn) ~ x), data = d)
for (m in c("epred", "predict")) {
  r <- tryCatch(fmt(suppressMessages(conditional_effects(
    fit, "x", method = m, int_conditions = list(x = xs)))[[1]]$estimate__),
    error = function(e) paste("ERROR", conditionMessage(e)))
  cat("cens(cn)", m, r, "\n")
}
