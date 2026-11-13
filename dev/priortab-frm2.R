# Messages, fits and round trips of lane priortab, before and after.
# PRIORTAB_LIB unset: rellib-r20 (before); set: the lane library.
# Rscript dev/priortab-frm2.R > dev/priortab-out/frm2-<arm>.txt
lib <- Sys.getenv("PRIORTAB_LIB")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
cat("frmtmb.sample", dirname(find.package("frmtmb.sample")), "\n")
options(width = 150, digits = 10)
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", conditionMessage(w), "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) {
           cat("MSG: ", conditionMessage(m))
           invokeRestart("muffleMessage")
         }),
       error = function(e) {
         cat("ERROR: ", conditionMessage(e), "\n")
         NULL
       })
  if (!is.null(r)) print(r)
  invisible(r)
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$y2 <- d$y + 0.5 * d$z

# messages
show("4312 list(sp1)", frm(y ~ z, data = d,
                           prior = list(set_prior("normal(0,1)"))))
show("4312 list(sp1, sp2)", frm(y ~ z, data = d,
  prior = list(set_prior("normal(0,1)"), set_prior("normal(0,1)",
                                                   class = "Intercept"))))
show("4312 list()", frm(y ~ z, data = d, prior = list()))
show("4312 list('flat')", frm(y ~ z, data = d, prior = list("flat")))
show("4316 prior = 'flat'", frm(y ~ z, data = d, prior = "flat"))
show("4316 prior = 'normal(0, 1)'", frm(y ~ z, data = d,
                                        prior = "normal(0, 1)"))
show("4283 one-column selection", {
  pe <- default_prior(y ~ z + x, data = d)
  pe$prior[pe$class == "b" & pe$coef == "z"] <- "normal(0, 1)"
  fixef(frm(y ~ z + x, data = d, prior = pe[, "prior", drop = FALSE]))
})
show("4283 two-column selection", {
  pe <- default_prior(y ~ z + x, data = d)
  fixef(frm(y ~ z + x, data = d, prior = pe[, c("prior", "class")]))
})
show("4283 row selection keeps working", {
  pe <- default_prior(y ~ z + x, data = d)
  pe$prior[pe$class == "b" & pe$coef == "z"] <- "normal(0, 1)"
  fixef(frm(y ~ z + x, data = d, prior = pe[pe$class == "b", ]))
})
show("4283 plain data frame with prior only (as.brmsprior)",
     as.brmsprior(data.frame(prior = "normal(0,1)", coef = "z")))

# 4925: uniform on a slope without bounds
show("4925 uniform(0, 0.2) no bounds", {
  f <- frm(y2 ~ z, data = d,
           prior = set_prior("uniform(0, 0.2)", class = "b", coef = "z"))
  c(fixef(f)[, 1], logLik = as.numeric(logLik(f)),
    conv = f$opt$convergence)
})
show("4925 uniform(0, 0.2) with bounds", {
  f <- frm(y2 ~ z, data = d,
           prior = set_prior("uniform(0, 0.2)", class = "b", coef = "z",
                             lb = 0, ub = 0.2))
  c(fixef(f)[, 1], logLik = as.numeric(logLik(f)),
    conv = f$opt$convergence)
})
show("4925 uniform(0, 2) interior, no bounds", {
  f <- frm(y2 ~ z, data = d,
           prior = set_prior("uniform(0, 2)", class = "b", coef = "z"))
  c(fixef(f)[, 1], logLik = as.numeric(logLik(f)),
    conv = f$opt$convergence)
})
show("4925 normal(0, 1), control", {
  f <- frm(y2 ~ z, data = d,
           prior = set_prior("normal(0, 1)", class = "b", coef = "z"))
  c(fixef(f)[, 1], logLik = as.numeric(logLik(f)),
    conv = f$opt$convergence)
})

# 4246: the frm_sample() announcement with a class b prior
f0 <- frm(y ~ z + x, data = d)
show("4246 announcement, class b prior", {
  invisible(frm_sample(f0, prior = set_prior("normal(0, 1)", class = "b"),
                       chains = 1, iter = 20, warmup = 10, seed = 1,
                       refresh = 0))
  "sampled"
})
show("4246 announcement, no prior", {
  invisible(frm_sample(f0, chains = 1, iter = 20, warmup = 10, seed = 1,
                       refresh = 0))
  "sampled"
})

# the passed-back tables: does any fit move?
show("round trip fit route (1 | g)", {
  pe <- default_prior(y ~ z + (1 | g), data = d)
  a <- frm(y ~ z + (1 | g), data = d)
  b <- frm(y ~ z + (1 | g), data = d, prior = pe)
  list(identical_par = identical(a$opt$par, b$opt$par),
       a = a$opt$objective, b = b$opt$objective,
       lower = b$lower, upper = b$upper)
})
show("round trip fit route student", {
  pe <- default_prior(y ~ z, data = d, family = student())
  a <- frm(y ~ z, data = d, family = student())
  b <- frm(y ~ z, data = d, family = student(), prior = pe)
  list(identical_par = identical(a$opt$par, b$opt$par),
       a = a$opt$objective, b = b$opt$objective)
})
show("round trip fit route zero_inflated_poisson", {
  dz <- transform(d, cnt = rpois(n, exp(0.5 + 0.3 * z)) * rbinom(n, 1, 0.7))
  pe <- default_prior(cnt ~ z, data = dz, family = zero_inflated_poisson())
  a <- frm(cnt ~ z, data = dz, family = zero_inflated_poisson())
  b <- frm(cnt ~ z, data = dz, family = zero_inflated_poisson(), prior = pe)
  list(identical_par = identical(a$opt$par, b$opt$par),
       a = a$opt$objective, b = b$opt$objective)
})
show("round trip sample route (1 | g): -log prior at the mode", {
  ps <- default_prior(y ~ z + (1 | g), data = d, route = "sample")
  fg <- frm(y ~ z + (1 | g), data = d)
  ri0 <- frmtmb:::resolve_prior_input(fg, frmtmb.sample:::default_priors_for(fg))
  ri1 <- frmtmb:::resolve_prior_input(fg, ps)
  p <- fg$estimates
  f0 <- frmtmb:::neg_log_prior_fn(ri0$entries)
  f1 <- frmtmb:::neg_log_prior_fn(ri1$entries)
  list(default = f0(p), table = f1(p), lower0 = ri0$lower, lower1 = ri1$lower)
})
show("round trip sample route student nu", {
  ps <- default_prior(y ~ z, data = d, family = student(), route = "sample")
  fs <- frm(y ~ z, data = d, family = student())
  ri0 <- frmtmb:::resolve_prior_input(fs, frmtmb.sample:::default_priors_for(fs))
  ri1 <- frmtmb:::resolve_prior_input(fs, ps)
  p <- fs$estimates
  list(default = frmtmb:::neg_log_prior_fn(ri0$entries)(p),
       table = frmtmb:::neg_log_prior_fn(ri1$entries)(p),
       lower0 = ri0$lower, lower1 = ri1$lower)
})
cat("\ndone\n")
