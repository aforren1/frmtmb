# Release review of 0.82.0, probe 3: NEWS bullets against the build, and
# gp_inner_pd() on an inner Hessian that is not finite.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
say <- function(label, expr) {
  conds <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    structure(list(msg = conditionMessage(e)), class = "probe_error")
  }), warning = function(w) {
    conds <<- c(conds, paste0("W[", class(w)[1], "]: ",
                              substr(conditionMessage(w), 1, 90)))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    conds <<- c(conds, paste0("M[", class(m)[1], "]: ",
                              substr(conditionMessage(m), 1, 90)))
    invokeRestart("muffleMessage")
  })
  cat("\n##", label, "\n")
  if (inherits(v, "probe_error")) cat("ERROR:", v$msg, "\n") else print(v)
  if (length(conds)) cat(paste0("  ", conds, "\n"), sep = "")
  invisible(v)
}

set.seed(11)
d <- data.frame(x = rnorm(120), z = rnorm(120),
                g = factor(rep(1:12, each = 10)))
d$y <- 1 + 0.5 * d$x + rnorm(12)[d$g] + 0.3 * d$x * rnorm(12)[d$g] +
  rnorm(120)

# NEWS: a fit that stops on a bound with an objective that is not finite
# is code 1 with the convergence warning (lognormal(0, 1), lb = 0, slope
# the data put below 0)
d$yn <- 1 - 2 * d$x + rnorm(120)
say("lognormal prior lb = 0 on a negative slope: code", {
  f <- frm(yn ~ x, data = d,
           prior = set_prior("lognormal(0, 1)", class = "b", coef = "x",
                             lb = 0))
  c(code = f$opt$convergence, objective = f$opt$objective,
    slope = unname(fixef(f)["x", "Estimate"]))
})

# NEWS: uniform(0, 0.2) on a slope with no lb/ub fits
say("uniform(0, 0.2) on a slope, no bounds", {
  f <- frm(y ~ x, data = d,
           prior = set_prior("uniform(0, 0.2)", class = "b", coef = "x"))
  c(code = f$opt$convergence, slope = unname(fixef(f)["x", "Estimate"]))
})

# NEWS: default_prior() + set_prior() combine
say("default_prior() + set_prior()", {
  p <- default_prior(y ~ x, data = d) + set_prior("normal(0, 1)", class = "b")
  c(class = class(p)[1], rows = nrow(p))
})

# NEWS: ranef() of a model without group-level terms stops
say("ranef() with no group terms", ranef(frm(y ~ x, data = d)))

# NEWS: ranef()$g holds every column; a block label returns that block
f2 <- frm(y ~ x + (1 | g) + (0 + x | g), data = d)
say("ranef(fit) names and $g columns", {
  r <- ranef(f2)
  list(names = names(r), g_cols = dimnames(r$g)[[3]],
       block = tryCatch(dimnames(r[["0 + x | g"]])[[3]],
                        error = function(e) conditionMessage(e)))
})
say("coef(fit)$g columns", colnames(coef(f2)$g[, , "Estimate"]))

# NEWS: predict(summary = FALSE) dimnames list(NULL, NULL)
say("predict(summary = FALSE) dimnames", {
  p <- predict(frm(y ~ x, data = d), summary = FALSE, ndraws = 5)
  dimnames(p)
})

# NEWS: Data: line for data passed by value
say("summary Data: line by do.call", {
  f <- do.call(frm, list(y ~ x, data = d))
  out <- capture.output(summary(f))
  grep("^ *Data:", out, value = TRUE)
})

# NEWS: df.residual() of y | mi() with 4 of 120 missing goes to 113
dm <- d
dm$y[c(3, 30, 60, 90)] <- NA
say("df.residual of y | mi() ~ x + z, 4 of 120 missing", {
  f <- frm(bf(y | mi() ~ x + z), data = dm)
  c(df_residual = df.residual(f), nobs = nobs(f),
    ll_nobs = attr(logLik(f), "nobs"))
})

# NEWS: the autocor gap message
dg <- data.frame(t = rep(1:6, 4), grp = factor(rep(1:4, each = 6)))
dg$y <- rnorm(24)
dg <- dg[!(dg$grp == "2" & dg$t == 3), ]
say("ar(cov = TRUE) with a skipped time level",
    logLik(frm(y ~ ar(time = t, gr = grp, cov = TRUE), data = dg)))

# NEWS: an ordinal mixture runs 10 starts by default, none with a start of
# your own, none in frm_allfit()
set.seed(3)
n <- 300
xo <- rnorm(n)
cls <- rbinom(n, 1, 0.5)
lat <- ifelse(cls == 1, 1.5 * xo, -1.5 * xo) + rlogis(n)
do <- data.frame(xo = xo, yo = cut(lat, c(-Inf, -1, 0, 1, Inf),
                                   labels = FALSE))
say("ordinal mixture default starts", {
  f <- frm(yo ~ xo, data = do, family = mixture(cumulative(), cumulative()))
  f$opt$mixture_starts
})

# gp_inner_pd() with a non-finite inner Hessian
library(Matrix)
fake <- list(fn = function(par) 0,
             env = list(random = 1:3, last.par = 0,
                        spHess = function(par, random) {
                          forceSymmetric(Matrix(diag(c(NaN, 2, 2)),
                                                sparse = TRUE))
                        }))
say("gp_inner_pd on a NaN inner Hessian (should be FALSE)",
    frmtmb:::gp_inner_pd(fake, 0, 1e-2))
fake$env$spHess <- function(par, random) {
  forceSymmetric(Matrix(diag(c(Inf, 2, 2)), sparse = TRUE))
}
say("gp_inner_pd on an Inf inner Hessian (should be FALSE)",
    frmtmb:::gp_inner_pd(fake, 0, 1e-2))
